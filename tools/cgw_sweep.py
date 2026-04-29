import json
import math
import sys
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed

import numpy as np
from scipy.linalg import expm


def vector_norm(v):
    return float(np.linalg.norm(v))


def add_vectors(a, b):
    return a + b


def subtract_vectors(a, b):
    return a - b


def scale_vector(v, scale):
    return v * scale


def dot(a, b):
    return float(np.dot(a, b))


def normalize_vector(v):
    norm = vector_norm(v)
    if not math.isfinite(norm) or norm < 1e-12:
        return v.copy()
    return v / norm


def solve_linear_system(matrix, rhs):
    try:
        return np.linalg.solve(np.array(matrix, dtype=float), np.array(rhs, dtype=float))
    except np.linalg.LinAlgError:
        return None


def determinant(matrix):
    try:
        return float(np.linalg.det(np.array(matrix, dtype=float)))
    except np.linalg.LinAlgError:
        return 0.0


def null_vector_square(matrix):
    a = np.array(matrix, dtype=float)
    _, _, vh = np.linalg.svd(a)
    return vh[-1, :]


def mat_mul(A, B):
    return np.matmul(np.array(A, dtype=float), np.array(B, dtype=float))


def mat_vec_mul(A, x):
    return np.matmul(np.array(A, dtype=float), np.array(x, dtype=float))


def transpose(A):
    return np.array(A, dtype=float).T


def invert_2x2(A):
    A = np.array(A, dtype=float)
    detA = A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0]
    if abs(detA) < 1e-12:
        return None
    return np.array([
        [A[1, 1] / detA, -A[0, 1] / detA],
        [-A[1, 0] / detA, A[0, 0] / detA],
    ])


def minimal_mass(q, params):
    delta = q[0] - q[1]
    coupling = params["m"] * params["b"] * math.cos(delta)
    return np.array([
        [params["m"] * params["b"] * params["b"], -coupling],
        [-coupling, params["M"] + params["m"] * (1 + params["a"] * params["a"])],
    ])


def gravity_vector(q, params):
    return np.array([
        params["m"] * params["b"] * math.sin(q[0]),
        -(params["M"] + params["m"] * (1 + params["a"])) * math.sin(q[1]),
    ])


def coriolis_vector(q, dq, params):
    term = params["m"] * params["b"] * math.sin(q[0] - q[1])
    return np.array([-term * dq[1] * dq[1], term * dq[0] * dq[0]])


def continuous_dynamics(x, params):
    q = np.array([x[0], x[1]])
    dq = np.array([x[2], x[3]])
    M = minimal_mass(q, params)
    c = coriolis_vector(q, dq, params)
    g = gravity_vector(q, params)
    rhs = np.array([-c[0] - g[0], -c[1] - g[1]])
    ddq = solve_linear_system(M, rhs)
    if ddq is None:
        return np.array([dq[0], dq[1], 0.0, 0.0])
    return np.array([dq[0], dq[1], ddq[0], ddq[1]])


def rk4_step(x, h, params):
    def add_scaled(base, delta, scale):
        return base + delta * scale

    k1 = continuous_dynamics(x, params)
    k2 = continuous_dynamics(add_scaled(x, k1, 0.5 * h), params)
    k3 = continuous_dynamics(add_scaled(x, k2, 0.5 * h), params)
    k4 = continuous_dynamics(add_scaled(x, k3, h), params)
    return x + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0


def integrate_passive(T, x0, params, steps):
    h = T / max(steps, 2)
    x = np.array(x0, dtype=float)
    trajectory = [x.copy()]
    times = [0.0]
    for i in range(steps):
        x = rk4_step(x, h, params)
        trajectory.append(x.copy())
        times.append((i + 1) * h)
    return {"xT": x, "trajectory": trajectory, "times": times}


def full_mass_matrix(sw_abs, st_abs, params):
    mb = params["m"] * params["b"]
    return np.array([
        [params["M"] + 2 * params["m"], 0, mb * math.cos(sw_abs), mb * math.cos(st_abs)],
        [0, params["M"] + 2 * params["m"], mb * math.sin(sw_abs), mb * math.sin(st_abs)],
        [mb * math.cos(sw_abs), mb * math.sin(sw_abs), params["m"] * params["b"] * params["b"], 0],
        [mb * math.cos(st_abs), mb * math.sin(st_abs), 0, params["m"] * params["b"] * params["b"]],
    ])


def impact_map(xT, gamma, params):
    sw_abs = xT[0] + gamma
    st_abs = xT[1] + gamma
    B = np.array([
        [0, -math.cos(st_abs)],
        [0, -math.sin(st_abs)],
        [1, 0],
        [0, 1],
    ])
    dq_full = mat_vec_mul(B, np.array([xT[2], xT[3]]))
    Mfb = full_mass_matrix(sw_abs, st_abs, params)
    inv_cols = []
    for e in np.eye(4):
        col = solve_linear_system(Mfb, e)
        if col is None:
            return np.array(xT, dtype=float)
        inv_cols.append(col)
    Minv = transpose(inv_cols)
    W = np.array([
        [1, 0],
        [0, 1],
        [math.cos(sw_abs), math.sin(sw_abs)],
        [0, 0],
    ])
    WT = transpose(W)
    Gd = mat_mul(mat_mul(WT, Minv), W)
    Gd_inv = invert_2x2(Gd)
    if Gd_inv is None:
        return np.array(xT, dtype=float)
    correction = mat_mul(mat_mul(Minv, W), mat_mul(Gd_inv, WT))
    I = np.eye(4)
    delta_core = I - correction
    swap = np.array([
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 0, 1],
        [0, 0, 1, 0],
    ])
    Delta = mat_mul(swap, delta_core)
    dq_plus = mat_vec_mul(Delta, dq_full)
    return np.array([xT[1], xT[0], dq_plus[2], dq_plus[3]])


def standstill_linearization(params):
    zero_state = np.zeros(4)
    A = numerical_jacobian(lambda x: continuous_dynamics(x, params), zero_state)
    impact_x = np.zeros((4, 4))
    gamma_impact = np.zeros(4)
    eps = 1e-6

    for j in range(4):
        xp = zero_state.copy()
        xm = zero_state.copy()
        xp[j] += eps
        xm[j] -= eps
        fp = impact_map(xp, 0.0, params)
        fm = impact_map(xm, 0.0, params)
        impact_x[:, j] = (fp - fm) / (2 * eps)

    gp = impact_map(zero_state, eps, params)
    gm = impact_map(zero_state, -eps, params)
    gamma_impact[:] = (gp - gm) / (2 * eps)

    return {
        "A": A,
        "impact_x": impact_x,
        "gamma_impact": gamma_impact,
    }


def linearized_standstill_jacobian(T, params, linearization):
    Phi = expm(linearization["A"] * T)
    Jx = linearization["impact_x"] @ Phi - np.eye(4)
    J = np.zeros((6, 6))
    J[0:4, 0:4] = Jx
    J[0:4, 4] = linearization["gamma_impact"]
    J[4, 0:4] = Phi[0, :] + Phi[1, :]
    J[4, 4] = 2.0
    J[5, 0:4] = 2.0 * Phi[0, :]
    J[5, 4] = 2.0
    J[5, 5] = -T
    return J


def res_passive(z, params):
    T = z[0]
    x0 = z[1:5]
    gamma = z[5]
    v_avg = z[6]
    if not (T > 0) or gamma < -0.05 or gamma > 0.28:
        return np.array([10, 10, 10, 10, 10, 10], dtype=float)
    flow = integrate_passive(T, x0, params, 120)
    xT = flow["xT"]
    xImpact = impact_map(xT, gamma, params)
    return np.array([
        xImpact[0] - x0[0],
        xImpact[1] - x0[1],
        xImpact[2] - x0[2],
        xImpact[3] - x0[3],
        xT[0] + xT[1] + 2 * gamma,
        2 * math.sin(xT[0] + gamma) - v_avg * T,
    ])


def numerical_jacobian(fn, x):
    base = fn(x)
    J = np.zeros((len(base), len(x)))
    for j in range(len(x)):
        eps = 1e-5 * max(1.0, abs(x[j]))
        xp = x.copy()
        xm = x.copy()
        xp[j] += eps
        xm[j] -= eps
        fp = fn(xp)
        fm = fn(xm)
        J[:, j] = (fp - fm) / (2 * eps)
    return J


def newton_solve(fn, x_init, max_iterations):
    x = np.array(x_init, dtype=float)
    for _ in range(max_iterations):
        r = fn(x)
        if vector_norm(r) < 1e-7:
            return {"x": x, "ok": True}
        J = numerical_jacobian(fn, x)
        delta = solve_linear_system(J, -r)
        if delta is None:
            return {"x": x, "ok": False}
        current_norm = vector_norm(r)
        improved = False
        step = 1.0
        while step > 1 / 64:
            candidate = add_vectors(x, scale_vector(delta, step))
            if vector_norm(fn(candidate)) < current_norm:
                x = candidate
                improved = True
                break
            step *= 0.5
        if not improved:
            return {"x": x, "ok": False}
    return {"x": x, "ok": vector_norm(fn(x)) < 1e-6}


def compute_continuation_tangent(z, params):
    J = numerical_jacobian(lambda inp: res_passive(inp, params), z)
    A = J[:, 0:6]
    b = -J[:, 6]
    u = solve_linear_system(A, b)
    if u is None:
        return None
    tangent = normalize_vector(np.concatenate([u, [1.0]]))
    return scale_vector(tangent, -1) if tangent[6] < 0 else tangent


def bifurcation_indicator(T, params, linearization):
    return determinant(linearized_standstill_jacobian(T, params, linearization))


def bisect_root(fn, left, right):
    a = left
    b = right
    fa = fn(a)
    for _ in range(24):
        mid = 0.5 * (a + b)
        fm = fn(mid)
        if abs(fm) < 1e-9:
            return mid
        if fa * fm <= 0:
            b = mid
        else:
            a = mid
            fa = fm
    return 0.5 * (a + b)


def find_bifurcation_periods(params, linearization):
    roots = []
    prev_T = 0.2
    prev_value = bifurcation_indicator(prev_T, params, linearization)
    T = 0.25
    while T <= 3.0 + 1e-12:
        value = bifurcation_indicator(T, params, linearization)
        if math.isfinite(prev_value) and math.isfinite(value) and prev_value * value < 0:
            roots.append(bisect_root(lambda tau: bifurcation_indicator(tau, params, linearization), prev_T, T))
        prev_T = T
        prev_value = value
        T += 0.1
    if len(roots) < 2:
        return [1.15, 2.05]
    return roots[:2]


def seed_from_bifurcation(Tbif, params, branch_index, linearization):
    base = np.array([Tbif, 0, 0, 0, 0, 0, 0], dtype=float)
    tangent_min = normalize_vector(null_vector_square(linearized_standstill_jacobian(Tbif, params, linearization)))
    if tangent_min[5] < 0:
        tangent_min = scale_vector(tangent_min, -1)
    step_scale = 0.035 if branch_index == 0 else 0.055
    seed = add_vectors(base, np.concatenate([[0.0], scale_vector(tangent_min, step_scale)]))
    fixed = newton_solve(lambda y: res_passive(np.array([y[0], y[1], y[2], y[3], y[4], y[5], seed[6]]), params), seed[:6], 8)
    if fixed["ok"]:
        return np.concatenate([fixed["x"], [seed[6]]])
    return np.array(
        [1.08, 0.11, -0.13, -0.24, 0.14, 0.012, 0.02] if branch_index == 0
        else [1.92, 0.22, -0.24, -0.42, 0.22, 0.018, 0.03],
        dtype=float,
    )


def solve_augmented_corrector(z_pred, tangent, params):
    return newton_solve(
        lambda z: np.concatenate([res_passive(z, params), [dot(tangent, subtract_vectors(z, z_pred))]]),
        z_pred,
        10,
    )


def compute_continuation_branch(seed, params, gamma_target_deg):
    z = np.array(seed, dtype=float)
    tangent = compute_continuation_tangent(z, params)
    if tangent is None:
        return []
    data = []
    count = 0
    gamma_target = gamma_target_deg * math.pi / 180.0

    def push_point(point):
        data.append({
            "T": point[0],
            "x0": point[1:5].copy(),
            "gamma": point[5],
            "gammaDeg": point[5] * 180.0 / math.pi,
            "vAvg": point[6],
        })

    push_point(z)
    while count < 80:
        if not math.isfinite(z[5]) or z[5] >= gamma_target:
            break
        count += 1
        z_pred = add_vectors(z, scale_vector(tangent, 0.03))
        corrected = solve_augmented_corrector(z_pred, tangent, params)
        if not corrected["ok"]:
            break
        if (not math.isfinite(corrected["x"][5]) or corrected["x"][5] < -1e-6):
            break
        if abs(corrected["x"][5] - z[5]) < 1e-6 and abs(corrected["x"][6] - z[6]) < 1e-6:
            break
        z = corrected["x"]
        tangent = compute_continuation_tangent(z, params)
        if tangent is None:
            break
        push_point(z)
    return data


def build_continuation_data(params, gamma_target_deg):
    linearization = standstill_linearization(params)
    bif = find_bifurcation_periods(params, linearization)
    return {
        "short": compute_continuation_branch(seed_from_bifurcation(bif[0], params, 0, linearization), params, gamma_target_deg),
        "long": compute_continuation_branch(seed_from_bifurcation(bif[1], params, 1, linearization), params, gamma_target_deg),
    }


def analyze_point(a, m, gamma_target_deg=15.0):
    params = {"a": a, "b": 1 - a, "m": m, "M": 1 - 2 * m}
    t0 = time.time()
    data = build_continuation_data(params, gamma_target_deg)
    elapsed = time.time() - t0
    out = {"a": a, "m": m, "elapsed_s": elapsed}
    for branch in ("short", "long"):
        series = data[branch]
        out[branch] = {
            "count": len(series),
            "max_gamma_deg": max((p["gammaDeg"] for p in series), default=float("nan")),
        }
    return out


def analyze_point_safe(a, m, gamma_target_deg=15.0):
    try:
        return analyze_point(a, m, gamma_target_deg)
    except Exception as exc:  # noqa: BLE001
        return {
            "a": a,
            "m": m,
            "elapsed_s": None,
            "error": f"{type(exc).__name__}: {exc}",
            "short": {"count": 0, "max_gamma_deg": float("nan")},
            "long": {"count": 0, "max_gamma_deg": float("nan")},
        }


def frange(start, stop, step):
    values = []
    current = start
    while current <= stop + 1e-12:
        values.append(round(current, 10))
        current += step
    return values


def sweep_grid(a_values, m_values, gamma_target_deg=15.0, workers=8):
    tasks = [(a, m, gamma_target_deg) for m in m_values for a in a_values]
    results = []
    executor_cls = ProcessPoolExecutor
    try:
        executor = executor_cls(max_workers=workers)
    except PermissionError:
        executor_cls = ThreadPoolExecutor
        executor = executor_cls(max_workers=workers)

    with executor:
        future_map = {
            executor.submit(analyze_point_safe, a, m, gamma_target_deg): (a, m)
            for a, m, gamma_target_deg in tasks
        }
        for future in as_completed(future_map):
            results.append(future.result())
    results.sort(key=lambda item: (item["m"], item["a"]))
    return results


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--grid":
        a_min = float(sys.argv[2])
        a_max = float(sys.argv[3])
        a_step = float(sys.argv[4])
        m_min = float(sys.argv[5])
        m_max = float(sys.argv[6])
        m_step = float(sys.argv[7])
        gamma = float(sys.argv[8]) if len(sys.argv) > 8 else 15.0
        workers = int(sys.argv[9]) if len(sys.argv) > 9 else 8
        data = sweep_grid(frange(a_min, a_max, a_step), frange(m_min, m_max, m_step), gamma, workers)
        print(json.dumps(data, indent=2))
    else:
        a = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
        m = float(sys.argv[2]) if len(sys.argv) > 2 else 0.25
        gamma = float(sys.argv[3]) if len(sys.argv) > 3 else 15.0
        print(json.dumps(analyze_point_safe(a, m, gamma), indent=2))

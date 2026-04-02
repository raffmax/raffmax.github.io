(function () {
  "use strict";

  function clamp(value, min, max) {
    return Math.min(max, Math.max(min, value));
  }

  function normalizeGammaDeg(value, min, max) {
    var clamped = clamp(value, min, max);
    return Math.abs(clamped) < 0.005 ? 0 : clamped;
  }

  function pointDistance(a, b) {
    return Math.hypot(a.x - b.x, a.y - b.y);
  }

  function distanceToSegment(point, start, end) {
    var dx = end.x - start.x;
    var dy = end.y - start.y;
    var lengthSquared = dx * dx + dy * dy;
    var t;
    var projection;
    if (lengthSquared < 1e-12) return pointDistance(point, start);
    t = ((point.x - start.x) * dx + (point.y - start.y) * dy) / lengthSquared;
    projection = {
      x: start.x + clamp(t, 0, 1) * dx,
      y: start.y + clamp(t, 0, 1) * dy
    };
    return pointDistance(point, projection);
  }

  function vectorNorm(v) {
    var sum = 0;
    for (var i = 0; i < v.length; i += 1) sum += v[i] * v[i];
    return Math.sqrt(sum);
  }

  function addVectors(a, b) {
    return a.map(function (value, index) { return value + b[index]; });
  }

  function subtractVectors(a, b) {
    return a.map(function (value, index) { return value - b[index]; });
  }

  function scaleVector(v, scale) {
    return v.map(function (value) { return value * scale; });
  }

  function dot(a, b) {
    var total = 0;
    for (var i = 0; i < a.length; i += 1) total += a[i] * b[i];
    return total;
  }

  function normalizeVector(v) {
    var norm = vectorNorm(v);
    if (!isFinite(norm) || norm < 1e-12) return v.slice();
    return v.map(function (value) { return value / norm; });
  }

  function solveLinearSystem(matrix, rhs) {
    var n = matrix.length;
    var a = [];
    var i;
    for (i = 0; i < n; i += 1) {
      a[i] = matrix[i].slice();
      a[i].push(rhs[i]);
    }

    for (var col = 0; col < n; col += 1) {
      var pivotRow = col;
      var pivotAbs = Math.abs(a[col][col]);
      for (var row = col + 1; row < n; row += 1) {
        var candidate = Math.abs(a[row][col]);
        if (candidate > pivotAbs) {
          pivotAbs = candidate;
          pivotRow = row;
        }
      }
      if (pivotAbs < 1e-10) return null;
      if (pivotRow !== col) {
        var tmp = a[col];
        a[col] = a[pivotRow];
        a[pivotRow] = tmp;
      }
      var pivot = a[col][col];
      for (var row2 = col + 1; row2 < n; row2 += 1) {
        var factor = a[row2][col] / pivot;
        for (var c = col; c <= n; c += 1) {
          a[row2][c] -= factor * a[col][c];
        }
      }
    }

    var x = new Array(n);
    for (var r = n - 1; r >= 0; r -= 1) {
      var sum = a[r][n];
      for (var c2 = r + 1; c2 < n; c2 += 1) {
        sum -= a[r][c2] * x[c2];
      }
      x[r] = sum / a[r][r];
    }
    return x;
  }

  function determinant(matrix) {
    var n = matrix.length;
    var a = matrix.map(function (row) { return row.slice(); });
    var sign = 1;
    var det = 1;
    for (var col = 0; col < n; col += 1) {
      var pivotRow = col;
      var pivotAbs = Math.abs(a[col][col]);
      for (var row = col + 1; row < n; row += 1) {
        var candidate = Math.abs(a[row][col]);
        if (candidate > pivotAbs) {
          pivotAbs = candidate;
          pivotRow = row;
        }
      }
      if (pivotAbs < 1e-12) return 0;
      if (pivotRow !== col) {
        var tmp = a[col];
        a[col] = a[pivotRow];
        a[pivotRow] = tmp;
        sign *= -1;
      }
      var pivot = a[col][col];
      det *= pivot;
      for (var row2 = col + 1; row2 < n; row2 += 1) {
        var factor = a[row2][col] / pivot;
        for (var c = col; c < n; c += 1) {
          a[row2][c] -= factor * a[col][c];
        }
      }
    }
    return sign * det;
  }

  function nullVectorSquare(matrix) {
    var rows = matrix.length;
    var cols = matrix[0].length;
    var a = matrix.map(function (row) { return row.slice(); });
    var pivotCols = [];
    var row = 0;

    for (var col = 0; col < cols && row < rows; col += 1) {
      var pivotRow = row;
      var pivotAbs = Math.abs(a[row][col]);
      for (var r = row + 1; r < rows; r += 1) {
        var candidate = Math.abs(a[r][col]);
        if (candidate > pivotAbs) {
          pivotAbs = candidate;
          pivotRow = r;
        }
      }
      if (pivotAbs < 1e-8) continue;
      if (pivotRow !== row) {
        var tmp = a[row];
        a[row] = a[pivotRow];
        a[pivotRow] = tmp;
      }
      var pivot = a[row][col];
      for (var c = col; c < cols; c += 1) {
        a[row][c] /= pivot;
      }
      for (var r2 = 0; r2 < rows; r2 += 1) {
        if (r2 === row) continue;
        var factor = a[r2][col];
        for (var c2 = col; c2 < cols; c2 += 1) {
          a[r2][c2] -= factor * a[row][c2];
        }
      }
      pivotCols.push(col);
      row += 1;
    }

    var freeCol = cols - 1;
    while (pivotCols.indexOf(freeCol) !== -1 && freeCol > 0) {
      freeCol -= 1;
    }

    var x = new Array(cols).fill(0);
    x[freeCol] = 1;
    for (var pr = pivotCols.length - 1; pr >= 0; pr -= 1) {
      var pc = pivotCols[pr];
      var sum = 0;
      for (var c3 = pc + 1; c3 < cols; c3 += 1) sum += a[pr][c3] * x[c3];
      x[pc] = -sum;
    }
    return x;
  }

  function matMul(A, B) {
    var rows = A.length;
    var cols = B[0].length;
    var inner = B.length;
    var out = [];
    for (var i = 0; i < rows; i += 1) {
      out[i] = [];
      for (var j = 0; j < cols; j += 1) {
        var sum = 0;
        for (var k = 0; k < inner; k += 1) sum += A[i][k] * B[k][j];
        out[i][j] = sum;
      }
    }
    return out;
  }

  function matAdd(A, B) {
    return A.map(function (row, i) {
      return row.map(function (value, j) {
        return value + B[i][j];
      });
    });
  }

  function matSub(A, B) {
    return A.map(function (row, i) {
      return row.map(function (value, j) {
        return value - B[i][j];
      });
    });
  }

  function matScale(A, scale) {
    return A.map(function (row) {
      return row.map(function (value) {
        return value * scale;
      });
    });
  }

  function identityMatrix(n) {
    var out = [];
    for (var i = 0; i < n; i += 1) {
      out[i] = [];
      for (var j = 0; j < n; j += 1) out[i][j] = i === j ? 1 : 0;
    }
    return out;
  }

  function matrix1Norm(A) {
    var cols = A[0].length;
    var maxCol = 0;
    for (var j = 0; j < cols; j += 1) {
      var sum = 0;
      for (var i = 0; i < A.length; i += 1) sum += Math.abs(A[i][j]);
      if (sum > maxCol) maxCol = sum;
    }
    return maxCol;
  }

  function matrixExponential(A) {
    var n = A.length;
    var norm = matrix1Norm(A);
    var squarings = norm > 0 ? Math.max(0, Math.ceil(Math.log(norm) / Math.log(2))) : 0;
    var scaled = matScale(A, 1 / Math.pow(2, squarings));
    var result = identityMatrix(n);
    var term = identityMatrix(n);

    for (var k = 1; k <= 24; k += 1) {
      term = matScale(matMul(term, scaled), 1 / k);
      result = matAdd(result, term);
      if (matrix1Norm(term) < 1e-12) break;
    }

    for (var s = 0; s < squarings; s += 1) result = matMul(result, result);
    return result;
  }

  function matVecMul(A, x) {
    return A.map(function (row) {
      var sum = 0;
      for (var i = 0; i < row.length; i += 1) sum += row[i] * x[i];
      return sum;
    });
  }

  function transpose(A) {
    var out = [];
    for (var i = 0; i < A[0].length; i += 1) {
      out[i] = [];
      for (var j = 0; j < A.length; j += 1) out[i][j] = A[j][i];
    }
    return out;
  }

  function flattenMatrix(A) {
    var out = [];
    for (var i = 0; i < A.length; i += 1) {
      for (var j = 0; j < A[i].length; j += 1) out.push(A[i][j]);
    }
    return out;
  }

  function unflattenMatrix(data, rows, cols, offset) {
    var out = [];
    var index = offset || 0;
    for (var i = 0; i < rows; i += 1) {
      out[i] = [];
      for (var j = 0; j < cols; j += 1) {
        out[i][j] = data[index];
        index += 1;
      }
    }
    return out;
  }

  function invert2x2(A) {
    var detA = A[0][0] * A[1][1] - A[0][1] * A[1][0];
    if (Math.abs(detA) < 1e-12) return null;
    return [
      [A[1][1] / detA, -A[0][1] / detA],
      [-A[1][0] / detA, A[0][0] / detA]
    ];
  }

  function makeSvgEl(tag) {
    return document.createElementNS("http://www.w3.org/2000/svg", tag);
  }

  function eventToElementCoordinates(event, element, width, height) {
    var rect = element.getBoundingClientRect();
    return {
      x: (event.clientX - rect.left) * width / Math.max(rect.width, 1),
      y: (event.clientY - rect.top) * height / Math.max(rect.height, 1)
    };
  }

  function normalizedPaperModel(params) {
    return {
      a: params.a,
      b: params.b,
      l0: params.a + params.b,
      ml: params.m,
      mh: params.M,
      g: 1
    };
  }

  // Eqs. (36)-(37): compass-gait dynamics in minimal coordinates q = [theta_sw, theta_st].
  function minimalMass(q, params) {
    var model = normalizedPaperModel(params);
    var alpha = q[1] - q[0];
    var coupling = model.ml * model.l0 * model.b * Math.cos(alpha);
    return [
      [model.ml * model.b * model.b, -coupling],
      [-coupling, (model.mh + model.ml) * model.l0 * model.l0 + model.ml * model.a * model.a]
    ];
  }

  function gravityVector(q, params) {
    var model = normalizedPaperModel(params);
    return [
      model.ml * model.b * model.g * Math.sin(q[0]),
      -(model.mh * model.l0 + model.ml * model.a + model.ml * model.l0) * model.g * Math.sin(q[1])
    ];
  }

  function coriolisMatrix(q, dq, params) {
    var model = normalizedPaperModel(params);
    var alpha = q[1] - q[0];
    var coupling = model.ml * model.l0 * model.b * Math.sin(alpha);
    return [
      [0, coupling * dq[1]],
      [-coupling * dq[0], 0]
    ];
  }

  function coriolisVector(q, dq, params) {
    return matVecMul(coriolisMatrix(q, dq, params), dq);
  }

  function continuousDynamics(x, params) {
    var q = [x[0], x[1]];
    var dq = [x[2], x[3]];
    var M = minimalMass(q, params);
    var c = coriolisVector(q, dq, params);
    var g = gravityVector(q, params);
    var rhs = [-c[0] - g[0], -c[1] - g[1]];
    var ddq = solveLinearSystem(M, rhs);
    if (!ddq) return [dq[0], dq[1], 0, 0];
    return [dq[0], dq[1], ddq[0], ddq[1]];
  }

  function continuousDynamicsJacobian(x, params) {
    var q = [x[0], x[1]];
    var dq = [x[2], x[3]];
    var model = normalizedPaperModel(params);
    var alpha = q[1] - q[0];
    var sinAlpha = Math.sin(alpha);
    var cosAlpha = Math.cos(alpha);
    var k = model.ml * model.l0 * model.b;
    var g1Scale = model.ml * model.b * model.g;
    var g2Scale = (model.mh * model.l0 + model.ml * model.a + model.ml * model.l0) * model.g;
    var M = minimalMass(q, params);
    var Minv = invert2x2(M);
    var s = k * sinAlpha;
    var c = k * cosAlpha;
    var N;
    var ddq;

    if (!Minv) {
      return [
        [0, 0, 1, 0],
        [0, 0, 0, 1],
        [0, 0, 0, 0],
        [0, 0, 0, 0]
      ];
    }

    N = [
      -s * dq[1] * dq[1] - g1Scale * Math.sin(q[0]),
      s * dq[0] * dq[0] + g2Scale * Math.sin(q[1])
    ];
    ddq = matVecMul(Minv, N);

    function columnFromDerivatives(dM, dN) {
      return matVecMul(Minv, subtractVectors(dN, matVecMul(dM, ddq)));
    }

    var dMdq1 = [
      [0, -s],
      [-s, 0]
    ];
    var dMdq2 = [
      [0, s],
      [s, 0]
    ];
    var dNdq1 = [
      c * dq[1] * dq[1] - g1Scale * Math.cos(q[0]),
      -c * dq[0] * dq[0]
    ];
    var dNdq2 = [
      -c * dq[1] * dq[1],
      c * dq[0] * dq[0] + g2Scale * Math.cos(q[1])
    ];
    var dNdv1 = [0, 2 * s * dq[0]];
    var dNdv2 = [-2 * s * dq[1], 0];
    var zeroMassDerivative = [
      [0, 0],
      [0, 0]
    ];
    var ddqQ1 = columnFromDerivatives(dMdq1, dNdq1);
    var ddqQ2 = columnFromDerivatives(dMdq2, dNdq2);
    var ddqV1 = columnFromDerivatives(zeroMassDerivative, dNdv1);
    var ddqV2 = columnFromDerivatives(zeroMassDerivative, dNdv2);

    return [
      [0, 0, 1, 0],
      [0, 0, 0, 1],
      [ddqQ1[0], ddqQ2[0], ddqV1[0], ddqV2[0]],
      [ddqQ1[1], ddqQ2[1], ddqV1[1], ddqV2[1]]
    ];
  }

  function rk4Step(x, h, params) {
    function addScaled(base, delta, scale) {
      return base.map(function (value, index) { return value + delta[index] * scale; });
    }
    var k1 = continuousDynamics(x, params);
    var k2 = continuousDynamics(addScaled(x, k1, 0.5 * h), params);
    var k3 = continuousDynamics(addScaled(x, k2, 0.5 * h), params);
    var k4 = continuousDynamics(addScaled(x, k3, h), params);
    return x.map(function (value, index) {
      return value + h * (k1[index] + 2 * k2[index] + 2 * k3[index] + k4[index]) / 6;
    });
  }

  function integratePassive(T, x0, params, steps) {
    var h = T / Math.max(steps, 2);
    var x = x0.slice();
    var trajectory = [x.slice()];
    var times = [0];
    for (var i = 0; i < steps; i += 1) {
      x = rk4Step(x, h, params);
      trajectory.push(x.slice());
      times.push((i + 1) * h);
    }
    return { xT: x, trajectory: trajectory, times: times };
  }

  function rk4StepAugmented(X, h, params) {
    function addScaled(base, delta, scale) {
      return base.map(function (value, index) { return value + delta[index] * scale; });
    }

    function augmentedDynamics(state) {
      var x = state.slice(0, 4);
      var phi = unflattenMatrix(state, 4, 4, 4);
      var f = continuousDynamics(x, params);
      var A = continuousDynamicsJacobian(x, params);
      return f.concat(flattenMatrix(matMul(A, phi)));
    }

    var k1 = augmentedDynamics(X);
    var k2 = augmentedDynamics(addScaled(X, k1, 0.5 * h));
    var k3 = augmentedDynamics(addScaled(X, k2, 0.5 * h));
    var k4 = augmentedDynamics(addScaled(X, k3, h));
    return X.map(function (value, index) {
      return value + h * (k1[index] + 2 * k2[index] + 2 * k3[index] + k4[index]) / 6;
    });
  }

  function integratePassiveWithSensitivity(T, x0, params, steps) {
    var h = T / Math.max(steps, 2);
    var X = x0.slice().concat(flattenMatrix(identityMatrix(4)));
    for (var i = 0; i < steps; i += 1) {
      X = rk4StepAugmented(X, h, params);
    }
    var xT = X.slice(0, 4);
    return {
      xT: xT,
      xT_T: continuousDynamics(xT, params),
      xT_x0: unflattenMatrix(X, 4, 4, 4)
    };
  }

  function impactVelocityMapData(xT, params) {
    var model = normalizedPaperModel(params);
    var alpha = xT[1] - xT[0];
    var cosAlpha = Math.cos(alpha);
    var sinAlpha = Math.sin(alpha);
    var mlab = model.ml * model.a * model.b;
    var h2 = model.mh * model.l0 * model.l0 + 2 * model.ml * model.a * model.l0;
    var l2 = model.ml * model.l0 * model.l0 + model.ml * model.a * model.a + model.mh * model.l0 * model.l0;
    var k = model.ml * model.b * model.l0;
    var mlbl0cos = k * cosAlpha;
    var qMinus = [
      [-mlab, h2 * cosAlpha - mlab],
      [0, -mlab]
    ];
    var qPlus = [
      [model.ml * model.b * model.b - mlbl0cos, l2 - mlbl0cos],
      [model.ml * model.b * model.b, -mlbl0cos]
    ];
    var qPlusInv = invert2x2(qPlus);
    if (!qPlusInv) return null;
    var velocityMap = matMul(qPlusInv, qMinus);
    var dQMinus = [
      [0, -h2 * sinAlpha],
      [0, 0]
    ];
    var dQPlus = [
      [k * sinAlpha, k * sinAlpha],
      [0, k * sinAlpha]
    ];
    return {
      velocityMap: velocityMap,
      dVelocityMap: matMul(qPlusInv, matSub(dQMinus, matMul(dQPlus, velocityMap)))
    };
  }

  // Eqs. (38)-(39): impact/reset map written directly in the appendix form.
  function impactMap(xT, gamma, params) {
    void gamma;
    var mapData = impactVelocityMapData(xT, params);
    if (!mapData) return xT.slice();
    var dqPlus = matVecMul(mapData.velocityMap, [xT[2], xT[3]]);
    return [xT[1], xT[0], dqPlus[0], dqPlus[1]];
  }

  function impactJacobian(xT, params) {
    var mapData = impactVelocityMapData(xT, params);
    var alphaSensitivity;
    if (!mapData) {
      return identityMatrix(4);
    }
    alphaSensitivity = matVecMul(mapData.dVelocityMap, [xT[2], xT[3]]);
    return [
      [0, 1, 0, 0],
      [1, 0, 0, 0],
      [-alphaSensitivity[0], alphaSensitivity[0], mapData.velocityMap[0][0], mapData.velocityMap[0][1]],
      [-alphaSensitivity[1], alphaSensitivity[1], mapData.velocityMap[1][0], mapData.velocityMap[1][1]]
    ];
  }

  function standstillLinearization(params) {
    var zeroState = [0, 0, 0, 0];
    return {
      A: continuousDynamicsJacobian(zeroState, params),
      impactX: impactJacobian(zeroState, params),
      gammaImpact: [0, 0, 0, 0]
    };
  }

  function linearizedStandstillJacobian(T, params, linearization) {
    var Phi = matrixExponential(matScale(linearization.A, T));
    var Jx = matSub(matMul(linearization.impactX, Phi), identityMatrix(4));
    var J = [];

    for (var i = 0; i < 4; i += 1) {
      J[i] = Jx[i].slice();
      J[i].push(linearization.gammaImpact[i]);
      J[i].push(0);
    }

    J[4] = [];
    J[5] = [];
    for (var j = 0; j < 4; j += 1) {
      J[4][j] = Phi[0][j] + Phi[1][j];
      J[5][j] = 2 * Phi[0][j];
    }
    J[4].push(2);
    J[4].push(0);
    J[5].push(2);
    J[5].push(-T);
    return J;
  }

  function resPassive(z, params) {
    var T = z[0];
    var x0 = z.slice(1, 5);
    var gamma = z[5];
    var vAvg = z[6];
    if (!(T > 0) || gamma < -0.05 || gamma > 0.28) return [10, 10, 10, 10, 10, 10];
    var flow = integratePassive(T, x0, params, 120);
    var xT = flow.xT;
    var xImpact = impactMap(xT, gamma, params);
    return [
      xImpact[0] - x0[0],
      xImpact[1] - x0[1],
      xImpact[2] - x0[2],
      xImpact[3] - x0[3],
      xT[0] + xT[1] + 2 * gamma,
      2 * Math.sin(xT[0] + gamma) - vAvg * T
    ];
  }

  function rowTimesMatrix(row, M) {
    return matVecMul(transpose(M), row);
  }

  function zeroMatrix(rows, cols) {
    var out = [];
    for (var i = 0; i < rows; i += 1) {
      out[i] = new Array(cols).fill(0);
    }
    return out;
  }

  function resPassiveWithJacobian(z, params) {
    var T = z[0];
    var x0 = z.slice(1, 5);
    var gamma = z[5];
    var vAvg = z[6];
    if (!(T > 0) || gamma < -0.05 || gamma > 0.28) {
      return {
        res: [10, 10, 10, 10, 10, 10],
        J: zeroMatrix(6, 7)
      };
    }

    var flow = integratePassiveWithSensitivity(T, x0, params, 120);
    var xT = flow.xT;
    var xImpact = impactMap(xT, gamma, params);
    var impactX = impactJacobian(xT, params);
    var topT = matVecMul(impactX, flow.xT_T);
    var topX0 = matSub(matMul(impactX, flow.xT_x0), identityMatrix(4));
    var phaseSensitivity = rowTimesMatrix([1, 1, 0, 0], flow.xT_x0);
    var speedFactor = 2 * Math.cos(xT[0] + gamma);
    var speedSensitivity = rowTimesMatrix([speedFactor, 0, 0, 0], flow.xT_x0);
    var J = [];

    for (var i = 0; i < 4; i += 1) {
      J.push([topT[i]].concat(topX0[i]).concat([0, 0]));
    }
    J.push([flow.xT_T[0] + flow.xT_T[1]].concat(phaseSensitivity).concat([2, 0]));
    J.push([speedFactor * flow.xT_T[0] - vAvg].concat(speedSensitivity).concat([speedFactor, -T]));

    return {
      res: [
        xImpact[0] - x0[0],
        xImpact[1] - x0[1],
        xImpact[2] - x0[2],
        xImpact[3] - x0[3],
        xT[0] + xT[1] + 2 * gamma,
        2 * Math.sin(xT[0] + gamma) - vAvg * T
      ],
      J: J
    };
  }

  function resFixedGamma(y, gamma, params) {
    var flow = integratePassive(y[0], y.slice(1, 5), params, 120);
    var xT = flow.xT;
    var xImpact = impactMap(xT, gamma, params);
    return [
      xImpact[0] - y[1],
      xImpact[1] - y[2],
      xImpact[2] - y[3],
      xImpact[3] - y[4],
      xT[0] + xT[1] + 2 * gamma
    ];
  }

  function resFixedGammaWithJacobian(y, gamma, params) {
    var T = y[0];
    var x0 = y.slice(1, 5);
    if (!(T > 0)) {
      return {
        res: [10, 10, 10, 10, 10],
        J: zeroMatrix(5, 5)
      };
    }

    var flow = integratePassiveWithSensitivity(T, x0, params, 120);
    var xT = flow.xT;
    var xImpact = impactMap(xT, gamma, params);
    var impactX = impactJacobian(xT, params);
    var topT = matVecMul(impactX, flow.xT_T);
    var topX0 = matSub(matMul(impactX, flow.xT_x0), identityMatrix(4));
    var phaseSensitivity = rowTimesMatrix([1, 1, 0, 0], flow.xT_x0);
    var J = [];

    for (var i = 0; i < 4; i += 1) {
      J.push([topT[i]].concat(topX0[i]));
    }
    J.push([flow.xT_T[0] + flow.xT_T[1]].concat(phaseSensitivity));

    return {
      res: [
        xImpact[0] - y[1],
        xImpact[1] - y[2],
        xImpact[2] - y[3],
        xImpact[3] - y[4],
        xT[0] + xT[1] + 2 * gamma
      ],
      J: J
    };
  }

  function newtonSolveWithJacobian(evaluate, xInit, maxIterations) {
    var x = xInit.slice();
    for (var iter = 0; iter < maxIterations; iter += 1) {
      var evaluation = evaluate(x);
      var r = evaluation.res;
      if (vectorNorm(r) < 1e-7) return { x: x, ok: true };
      var delta = solveLinearSystem(evaluation.J, r.map(function (value) { return -value; }));
      if (!delta) return { x: x, ok: false };
      var currentNorm = vectorNorm(r);
      var improved = false;
      for (var step = 1; step > 1 / 64; step *= 0.5) {
        var candidate = addVectors(x, scaleVector(delta, step));
        if (vectorNorm(evaluate(candidate).res) < currentNorm) {
          x = candidate;
          improved = true;
          break;
        }
      }
      if (!improved) return { x: x, ok: false };
    }
    return { x: x, ok: vectorNorm(evaluate(x).res) < 1e-6 };
  }

  function computeContinuationTangent(z, params) {
    var J = resPassiveWithJacobian(z, params).J;
    var A = J.map(function (row) { return row.slice(0, 6); });
    var b = J.map(function (row) { return -row[6]; });
    var u = solveLinearSystem(A, b);
    if (!u) return null;
    var tangent = normalizeVector(u.concat([1]));
    return tangent[6] < 0 ? scaleVector(tangent, -1) : tangent;
  }

  function bifurcationIndicator(T, params, linearization) {
    return determinant(linearizedStandstillJacobian(T, params, linearization));
  }

  function bisectRoot(fn, left, right) {
    var a = left;
    var b = right;
    var fa = fn(a);
    for (var i = 0; i < 24; i += 1) {
      var mid = 0.5 * (a + b);
      var fm = fn(mid);
      if (Math.abs(fm) < 1e-9) return mid;
      if (fa * fm <= 0) {
        b = mid;
      } else {
        a = mid;
        fa = fm;
      }
    }
    return 0.5 * (a + b);
  }

  function findBifurcationPeriods(params, linearization) {
    var roots = [];
    var prevT = 0.2;
    var prevValue = bifurcationIndicator(prevT, params, linearization);
    for (var T = 0.25; T <= 3.0; T += 0.1) {
      var value = bifurcationIndicator(T, params, linearization);
      if (isFinite(prevValue) && isFinite(value) && prevValue * value < 0) {
        roots.push(bisectRoot(function (tau) {
          return bifurcationIndicator(tau, params, linearization);
        }, prevT, T));
      }
      prevT = T;
      prevValue = value;
    }
    if (roots.length < 2) return [1.15, 2.05];
    return roots.slice(0, 2);
  }

  function seedFromBifurcation(Tbif, params, branchIndex, linearization) {
    var base = [Tbif, 0, 0, 0, 0, 0, 0];
    var tangentMin = normalizeVector(nullVectorSquare(linearizedStandstillJacobian(Tbif, params, linearization)));
    if (tangentMin[5] < 0) tangentMin = scaleVector(tangentMin, -1);
    var scale = branchIndex === 0 ? 0.035 : 0.055;
    var seed = addVectors(base, [0].concat(scaleVector(tangentMin, scale)));
    var fixed = newtonSolveWithJacobian(function (y) {
      var evaluation = resPassiveWithJacobian([y[0], y[1], y[2], y[3], y[4], y[5], seed[6]], params);
      return {
        res: evaluation.res,
        J: evaluation.J.map(function (row) { return row.slice(0, 6); })
      };
    }, seed.slice(0, 6), 8);
    if (fixed.ok) return fixed.x.concat([seed[6]]);
    return branchIndex === 0
      ? [1.08, 0.11, -0.13, -0.24, 0.14, 0.012, 0.02]
      : [1.92, 0.22, -0.24, -0.42, 0.22, 0.018, 0.03];
  }

  function solveAugmentedCorrector(zPred, tangent, params) {
    return newtonSolveWithJacobian(function (z) {
      var evaluation = resPassiveWithJacobian(z, params);
      evaluation.res.push(dot(tangent, subtractVectors(z, zPred)));
      evaluation.J.push(tangent.slice());
      return evaluation;
    }, zPred, 10);
  }

  function computeContinuationBranch(seed, params, gammaTargetDeg) {
    var z = seed.slice();
    var tangent = computeContinuationTangent(z, params);
    if (!tangent) return [];
    var data = [];
    var count = 0;
    var gammaTarget = gammaTargetDeg * Math.PI / 180;

    function pushPoint(point) {
      data.push({
        T: point[0],
        x0: point.slice(1, 5),
        gamma: point[5],
        gammaDeg: point[5] * 180 / Math.PI,
        vAvg: point[6]
      });
    }

    pushPoint(z);
    while (count < 80) {
      if (!isFinite(z[5])) break;
      if (z[5] >= gammaTarget) break;
      count += 1;
      var zPred = addVectors(z, scaleVector(tangent, 0.03));
      var corrected = solveAugmentedCorrector(zPred, tangent, params);
      if (!corrected.ok) break;
      if (!isFinite(corrected.x[5]) || corrected.x[5] < -1e-6) break;
      if (Math.abs(corrected.x[5] - z[5]) < 1e-6 && Math.abs(corrected.x[6] - z[6]) < 1e-6) break;
      z = corrected.x;
      tangent = computeContinuationTangent(z, params);
      if (!tangent) break;
      pushPoint(z);
    }
    return data;
  }

  function buildContinuationData(params, gammaTargetDeg) {
    var linearization = standstillLinearization(params);
    var bif = findBifurcationPeriods(params, linearization);
    return {
      short: computeContinuationBranch(seedFromBifurcation(bif[0], params, 0, linearization), params, gammaTargetDeg),
      long: computeContinuationBranch(seedFromBifurcation(bif[1], params, 1, linearization), params, gammaTargetDeg)
    };
  }

  function closestBranchPoint(branchData, gammaDeg) {
    if (!branchData || !branchData.length) return null;
    var best = branchData[0];
    var bestDistance = Math.abs(best.gammaDeg - gammaDeg);
    for (var i = 1; i < branchData.length; i += 1) {
      var distance = Math.abs(branchData[i].gammaDeg - gammaDeg);
      if (distance < bestDistance) {
        best = branchData[i];
        bestDistance = distance;
      }
    }
    return best;
  }

  function solveOrbitForGamma(gammaDeg, branchName, params, continuation) {
    var branchData = continuation && continuation[branchName];
    if (!branchData || !branchData.length) return null;
    if (Math.abs(gammaDeg) < 1e-12) {
      return {
        T: 1,
        x0: [0, 0, 0, 0],
        xT: [0, 0, 0, 0],
        xImpact: [0, 0, 0, 0],
        gamma: 0,
        gammaDeg: 0,
        solveGamma: 0,
        solveGammaDeg: 0,
        vAvg: 0,
        branch: branchName,
        trajectory: [[0, 0, 0, 0], [0, 0, 0, 0]],
        times: [0, 1],
        isStandstill: true
      };
    }
    var gammaSolveDeg = Math.abs(gammaDeg);
    var gammaSolve = gammaSolveDeg * Math.PI / 180;
    var gammaDisplay = gammaDeg * Math.PI / 180;
    var seed = closestBranchPoint(branchData, gammaSolveDeg);
    var yInit = [seed.T].concat(seed.x0);
    var solved = newtonSolveWithJacobian(function (y) {
      return resFixedGammaWithJacobian(y, gammaSolve, params);
    }, yInit, 8);
    var y = solved.ok ? solved.x : yInit;
    var flow = integratePassive(y[0], y.slice(1, 5), params, 160);
    var xT = flow.xT;
    var xImpact = impactMap(xT, gammaSolve, params);
    return {
      T: y[0],
      x0: y.slice(1, 5),
      xT: xT,
      xImpact: xImpact,
      gamma: gammaDisplay,
      gammaDeg: gammaDeg,
      solveGamma: gammaSolve,
      solveGammaDeg: gammaSolveDeg,
      vAvg: 2 * Math.sin(xT[0] + gammaSolve) / y[0],
      branch: branchName,
      trajectory: flow.trajectory,
      times: flow.times,
      isStandstill: false
    };
  }

  document.addEventListener("DOMContentLoaded", function () {
    var canvas = document.getElementById("cgw-canvas");
    var bifurcationSvg = document.getElementById("cgw-bifurcation");
    var interactive = document.getElementById("compass-gait-explorer");
    var loadingOverlay = document.getElementById("cgw-loading");
    var loadingMessage = document.getElementById("cgw-loading-message");
    if (!canvas || !bifurcationSvg || !interactive || !loadingOverlay || !loadingMessage) return;

    var ctx = canvas.getContext("2d");
    var bifurcationSize = { width: 820, height: 460 };
    var gammaRange = { min: -15, max: 15 };
    var defaultCanvasFontFamily = (window.getComputedStyle(document.body).fontFamily || "-apple-system, \"Segoe UI\", Arial, sans-serif").trim();
    var controls = {
      a: document.getElementById("cgw-a"),
      m: document.getElementById("cgw-m"),
      aValue: document.getElementById("cgw-a-value"),
      mValue: document.getElementById("cgw-m-value"),
      gammaValue: document.getElementById("cgw-gamma-value"),
      branch: Array.prototype.slice.call(document.querySelectorAll("input[name='cgw-branch']"))
    };

    var matlabColors = {
      stance: "rgb(0,65,145)",
      swing: "rgb(50,50,50)",
      mass: "rgb(0,190,255)",
      ground: "#1d1d1d",
      short: "#2f6f3e",
      long: "#7b4aa1"
    };

    var state = {
      a: parseFloat(controls.a.value),
      m: parseFloat(controls.m.value),
      gammaDeg: 2,
      branch: "short",
      time: 0,
      lastFrame: null,
      continuation: null,
      currentOrbit: null,
      paramKey: "",
      loadingToken: 0,
      groundMarkerOffsetScreen: 0,
      groundMarkerLastAnchorWorldX: null,
      bifurcationHitPoints: [],
      slopeInteraction: null,
      slopeDrag: null,
      orbitUpdateScheduled: false,
      lastCanvasPointer: null,
      lastBifurcationPointer: null
    };

    function currentParams() {
      return {
        a: state.a,
        b: 1 - state.a,
        m: state.m,
        M: 1 - 2 * state.m
      };
    }

    function syncStateFromControls() {
      state.a = parseFloat(controls.a.value);
      state.m = parseFloat(controls.m.value);
      controls.branch.forEach(function (input) {
        if (input.checked) state.branch = input.value;
      });
      state.gammaDeg = normalizeGammaDeg(state.gammaDeg, gammaRange.min, gammaRange.max);
    }

    function renderControlValues() {
      controls.aValue.innerHTML = state.a.toFixed(2) + " <i>l</i><sub>0</sub>";
      controls.mValue.innerHTML = state.m.toFixed(2) + " <i>m</i><sub>0</sub>";
      state.gammaDeg = normalizeGammaDeg(state.gammaDeg, gammaRange.min, gammaRange.max);
      controls.gammaValue.innerHTML = state.gammaDeg.toFixed(2) + "&deg;";
      controls.branch.forEach(function (input) {
        input.checked = input.value === state.branch;
      });
    }

    function setLoading(isLoading, message) {
      interactive.classList.toggle("is-loading", isLoading);
      loadingOverlay.setAttribute("aria-hidden", isLoading ? "false" : "true");
      if (message) loadingMessage.textContent = message;
    }

    function gammaAbsMax() {
      return Math.max(Math.abs(gammaRange.min), Math.abs(gammaRange.max));
    }

    function rebuildContinuationIfNeeded() {
      var gammaMax = gammaAbsMax();
      var key = state.a.toFixed(3) + "|" + state.m.toFixed(3) + "|" + gammaMax.toFixed(2);
      if (state.continuation && state.paramKey === key) return;
      state.paramKey = key;
      state.continuation = buildContinuationData(currentParams(), gammaMax);
    }

    function worldToScreen(x, y, viewport) {
      return {
        x: viewport.originX + x * viewport.scale,
        y: viewport.originY - y * viewport.scale
      };
    }

    function drawCircle(point, radius, fillStyle) {
      ctx.beginPath();
      ctx.arc(point.x, point.y, radius, 0, 2 * Math.PI);
      ctx.fillStyle = fillStyle;
      ctx.fill();
    }

    function lerpState(a, b, alpha) {
      return a.map(function (value, index) {
        return value + (b[index] - value) * alpha;
      });
    }

    function drawWalker() {
      var width = canvas.width;
      var height = canvas.height;
      ctx.clearRect(0, 0, width, height);
      if (!state.currentOrbit) {
        state.slopeInteraction = null;
        return;
      }

      var orbit = state.currentOrbit;
      var params = currentParams();
      var phaseTime = (((state.time + 0.18 * orbit.T) % orbit.T) + orbit.T) % orbit.T;
      var samplePosition = phaseTime / orbit.T * (orbit.trajectory.length - 1);
      var frameIndex = Math.floor(samplePosition);
      var frameAlpha = samplePosition - frameIndex;
      var xCurrent = orbit.trajectory[Math.min(frameIndex, orbit.trajectory.length - 1)];
      var xNext = frameIndex + 1 < orbit.trajectory.length ? orbit.trajectory[frameIndex + 1] : orbit.xImpact;
      var x = lerpState(xCurrent, xNext, frameAlpha);
      var walkDirection = orbit.gammaDeg < 0 ? -1 : 1;
      var stride = 2 * Math.sin(orbit.xT[0] + orbit.solveGamma);
      var stepIndex = Math.floor(state.time / orbit.T);
      var stanceX = walkDirection * stepIndex * stride;
      var groundSlope = -Math.tan(orbit.gamma);

      function groundY(xPos) {
        return groundSlope * xPos;
      }

      var stanceFoot = { x: stanceX, y: groundY(stanceX) };
      var hip = {
        x: stanceX - walkDirection * Math.sin(x[1]),
        y: groundY(stanceX) + Math.cos(x[1])
      };
      var swingFoot = {
        x: hip.x + walkDirection * Math.sin(x[0]),
        y: hip.y - Math.cos(x[0])
      };
      var stanceCom = {
        x: stanceFoot.x + params.a * (hip.x - stanceFoot.x),
        y: stanceFoot.y + params.a * (hip.y - stanceFoot.y)
      };
      var swingCom = {
        x: swingFoot.x + params.a * (hip.x - swingFoot.x),
        y: swingFoot.y + params.a * (hip.y - swingFoot.y)
      };
      var footMidpoint = {
        x: 0.5 * (stanceFoot.x + swingFoot.x),
        y: 0.5 * (stanceFoot.y + swingFoot.y)
      };
      var scale = 450;
      var walkerAnchor = {
        x: width * 0.5,
        y: height * 0.82
      };
      var visualScale = scale / 250;
      var branchColor = state.branch === "long" ? matlabColors.long : matlabColors.short;

      ctx.fillStyle = "#ffffff";
      ctx.fillRect(0, 0, width, height);
      var groundAnchorScreen = {
        x: walkerAnchor.x,
        y: walkerAnchor.y
      };
      var groundAnchorWorldX = footMidpoint.x;
      var groundAnchorWorldY = groundY(groundAnchorWorldX);

      function groundToScreen(xPos) {
        return {
          x: groundAnchorScreen.x + (xPos - groundAnchorWorldX) * scale,
          y: groundAnchorScreen.y - (groundY(xPos) - groundAnchorWorldY) * scale
        };
      }

      var stanceFootScreen = groundToScreen(stanceFoot.x);

      function walkerToScreen(point) {
        return {
          x: stanceFootScreen.x + (point.x - stanceFoot.x) * scale,
          y: stanceFootScreen.y - (point.y - stanceFoot.y) * scale
        };
      }

      var swingFootScreen = walkerToScreen(swingFoot);
      var hipScreen = walkerToScreen(hip);
      var stanceComScreen = walkerToScreen(stanceCom);
      var swingComScreen = walkerToScreen(swingCom);

      ctx.strokeStyle = branchColor;
      ctx.lineWidth = 2.6 * visualScale;
      ctx.beginPath();
      var leftGround = groundToScreen(groundAnchorWorldX - 0.7);
      var rightGround = groundToScreen(groundAnchorWorldX + 0.7);
      ctx.moveTo(leftGround.x, leftGround.y);
      ctx.lineTo(rightGround.x, rightGround.y);
      ctx.stroke();

      ctx.strokeStyle = state.branch === "long" ? "rgba(123,74,161,0.28)" : "rgba(47,111,62,0.28)";
      ctx.lineWidth = 1.5 * visualScale;
      var tickSpacingWorld = 0.45;
      var tickSpacingScreen = tickSpacingWorld * scale;
      var tickHeight = 12 * visualScale;
      var groundScreenSlope = (rightGround.y - leftGround.y) / Math.max(rightGround.x - leftGround.x, 1e-9);

      function groundScreenYAtX(screenX) {
        return leftGround.y + (screenX - leftGround.x) * groundScreenSlope;
      }

      ctx.lineCap = "butt";
      if (state.groundMarkerLastAnchorWorldX == null) {
        state.groundMarkerLastAnchorWorldX = groundAnchorWorldX;
        state.groundMarkerOffsetScreen = 0;
      } else {
        state.groundMarkerOffsetScreen -= (groundAnchorWorldX - state.groundMarkerLastAnchorWorldX) * scale;
        state.groundMarkerLastAnchorWorldX = groundAnchorWorldX;
      }
      state.groundMarkerOffsetScreen = ((state.groundMarkerOffsetScreen % tickSpacingScreen) + tickSpacingScreen) % tickSpacingScreen;

      for (var tickScreenX = leftGround.x + state.groundMarkerOffsetScreen; tickScreenX <= rightGround.x; tickScreenX += tickSpacingScreen) {
        var snappedX = Math.round(tickScreenX) + 0.5;
        if (snappedX < leftGround.x || snappedX > rightGround.x) continue;
        var baseY = groundScreenYAtX(snappedX);
        ctx.beginPath();
        ctx.moveTo(snappedX, baseY - tickHeight);
        ctx.lineTo(snappedX, baseY);
        ctx.stroke();
      }

      var slopeHandleRadius = 9 * visualScale;
      function drawSlopeHandle(point) {
        ctx.fillStyle = "rgba(255,255,255,0.96)";
        ctx.strokeStyle = branchColor;
        ctx.lineWidth = 2 * visualScale;
        ctx.beginPath();
        ctx.arc(point.x, point.y, slopeHandleRadius, 0, 2 * Math.PI);
        ctx.fill();
        ctx.stroke();
      }

      drawSlopeHandle(leftGround);
      drawSlopeHandle(rightGround);

      state.slopeInteraction = {
        midpoint: groundAnchorScreen,
        leftHandle: leftGround,
        rightHandle: rightGround,
        handleRadius: slopeHandleRadius,
        lineRadius: 14 * visualScale
      };

      ctx.fillStyle = "#4d5b4f";
      ctx.font = (14 * visualScale).toFixed(1) + "px " + defaultCanvasFontFamily;

      ctx.lineCap = "round";
      ctx.strokeStyle = matlabColors.stance;
      ctx.lineWidth = 10 * visualScale;
      ctx.beginPath();
      ctx.moveTo(stanceFootScreen.x, stanceFootScreen.y);
      ctx.lineTo(hipScreen.x, hipScreen.y);
      ctx.stroke();

      ctx.strokeStyle = matlabColors.swing;
      ctx.beginPath();
      ctx.moveTo(hipScreen.x, hipScreen.y);
      ctx.lineTo(swingFootScreen.x, swingFootScreen.y);
      ctx.stroke();

      drawCircle(stanceComScreen, (7 + 24 * params.m) * visualScale, matlabColors.mass);
      drawCircle(swingComScreen, (7 + 24 * params.m) * visualScale, matlabColors.mass);
      drawCircle(hipScreen, (12 + 36 * params.M) * visualScale, matlabColors.mass);

      var stanceLegDx = hipScreen.x - stanceFootScreen.x;
      var stanceLegDy = hipScreen.y - stanceFootScreen.y;
      var stanceLegNorm = Math.max(Math.hypot(stanceLegDx, stanceLegDy), 1e-9);
      var stanceLegTangent = {
        x: stanceLegDx / stanceLegNorm,
        y: stanceLegDy / stanceLegNorm
      };
      var stanceLegNormal = {
        x: -stanceLegTangent.y,
        y: stanceLegTangent.x
      };
      if (stanceLegNormal.x < 0) {
        stanceLegNormal.x *= -1;
        stanceLegNormal.y *= -1;
      }
      var abOffset = 22 * visualScale;
      var aLabelPosition = {
        x: (stanceFootScreen.x + stanceComScreen.x) * 0.5 + stanceLegNormal.x * abOffset + stanceLegTangent.x * 2 * visualScale,
        y: (stanceFootScreen.y + stanceComScreen.y) * 0.5 + stanceLegNormal.y * abOffset + stanceLegTangent.y * 2 * visualScale
      };
      var bLabelPosition = {
        x: (hipScreen.x + stanceComScreen.x) * 0.5 + stanceLegNormal.x * abOffset - stanceLegTangent.x * 2 * visualScale,
        y: (hipScreen.y + stanceComScreen.y) * 0.5 + stanceLegNormal.y * abOffset - stanceLegTangent.y * 2 * visualScale
      };

      ctx.fillStyle = "#243027";
      ctx.font = (21 * visualScale).toFixed(1) + "px " + defaultCanvasFontFamily;
      ctx.save();
      ctx.textAlign = "center";
      ctx.textBaseline = "middle";
      ctx.fillText("M", hipScreen.x, hipScreen.y);
      ctx.fillText("m", stanceComScreen.x, stanceComScreen.y);
      ctx.fillText("m", swingComScreen.x, swingComScreen.y);
      ctx.textAlign = "center";
      ctx.textBaseline = "middle";
      ctx.fillText("a", aLabelPosition.x, aLabelPosition.y);
      ctx.fillText("b", bLabelPosition.x, bLabelPosition.y);
      ctx.restore();

      ctx.fillStyle = "#465247";
      ctx.font = "15px " + defaultCanvasFontFamily;
      updateCanvasCursor(state.lastCanvasPointer);
    }

    function drawBifurcationDiagram() {
      while (bifurcationSvg.firstChild) bifurcationSvg.removeChild(bifurcationSvg.firstChild);
      state.bifurcationHitPoints = [];
      var width = bifurcationSize.width;
      var height = bifurcationSize.height;
      var margin = { top: 28, right: 20, bottom: 68, left: 118 };
      var innerWidth = width - margin.left - margin.right;
      var innerHeight = height - margin.top - margin.bottom;
      var gammaMax = gammaAbsMax();
      var allPoints = [];
      if (state.continuation) {
        if (state.continuation.short) allPoints = allPoints.concat(state.continuation.short);
        if (state.continuation.long) allPoints = allPoints.concat(state.continuation.long);
      }
      if (state.currentOrbit && !state.currentOrbit.isStandstill) {
        allPoints.push({ T: state.currentOrbit.T });
      }
      var tMin = 0.9;
      var tMax = 2.1;
      if (allPoints.length) {
        tMin = allPoints[0].T;
        tMax = allPoints[0].T;
        allPoints.forEach(function (point) {
          tMin = Math.min(tMin, point.T);
          tMax = Math.max(tMax, point.T);
        });
        var pad = Math.max(0.08, 0.08 * (tMax - tMin || 1));
        tMin = Math.max(0, tMin - pad);
        tMax = tMax + pad;
      }
      var yTicks = 5;

      function xScale(value) {
        return margin.left + (value + gammaMax) / (2 * gammaMax) * innerWidth;
      }

      function yScale(value) {
        if (tMax === tMin) return margin.top + innerHeight * 0.5;
        return margin.top + innerHeight - (value - tMin) / (tMax - tMin) * innerHeight;
      }

      var bg = makeSvgEl("rect");
      bg.setAttribute("x", "0");
      bg.setAttribute("y", "0");
      bg.setAttribute("width", width);
      bg.setAttribute("height", height);
      bg.setAttribute("fill", "#fbfdf9");
      bifurcationSvg.appendChild(bg);

      for (var gx = -gammaMax; gx <= gammaMax; gx += 5) {
        var xTick = makeSvgEl("text");
        xTick.setAttribute("x", xScale(gx));
        xTick.setAttribute("y", margin.top + innerHeight + 35);
        xTick.setAttribute("text-anchor", "middle");
        xTick.setAttribute("fill", "#4e5c51");
        xTick.setAttribute("font-size", "28");
        xTick.textContent = gx.toFixed(0);
        bifurcationSvg.appendChild(xTick);

        var xTickMark = makeSvgEl("line");
        xTickMark.setAttribute("x1", xScale(gx));
        xTickMark.setAttribute("x2", xScale(gx));
        xTickMark.setAttribute("y1", margin.top + innerHeight);
        xTickMark.setAttribute("y2", margin.top + innerHeight + 7);
        xTickMark.setAttribute("stroke", "#425044");
        xTickMark.setAttribute("stroke-width", "1.5");
        bifurcationSvg.appendChild(xTickMark);
      }

      var standstillLabelY = margin.top + 16;
      var standstillLine = makeSvgEl("line");
      standstillLine.setAttribute("x1", xScale(0));
      standstillLine.setAttribute("x2", xScale(0));
      standstillLine.setAttribute("y1", standstillLabelY + 10);
      standstillLine.setAttribute("y2", margin.top + innerHeight);
      standstillLine.setAttribute("stroke", "#6f7b70");
      standstillLine.setAttribute("stroke-width", "2.5");
      standstillLine.setAttribute("stroke-dasharray", "7 6");
      standstillLine.style.cursor = "pointer";
      standstillLine.addEventListener("click", function () {
        state.gammaDeg = 0;
        updateOrbitOnly();
      });
      bifurcationSvg.appendChild(standstillLine);

      var standstillLabel = makeSvgEl("text");
      standstillLabel.setAttribute("x", xScale(0));
      standstillLabel.setAttribute("y", standstillLabelY);
      standstillLabel.setAttribute("text-anchor", "middle");
      standstillLabel.setAttribute("fill", "#5f6b60");
      standstillLabel.setAttribute("font-size", "28");
      standstillLabel.setAttribute("font-weight", "700");
      standstillLabel.textContent = "stand-still";
      standstillLabel.style.cursor = "pointer";
      standstillLabel.addEventListener("click", function () {
        state.gammaDeg = 0;
        updateOrbitOnly();
      });
      bifurcationSvg.appendChild(standstillLabel);

      for (var gyIndex = 0; gyIndex <= yTicks; gyIndex += 1) {
        var gy = tMin + (gyIndex / yTicks) * (tMax - tMin);
        var yTick = makeSvgEl("text");
        yTick.setAttribute("x", margin.left - 15);
        yTick.setAttribute("y", yScale(gy) + 5);
        yTick.setAttribute("text-anchor", "end");
        yTick.setAttribute("fill", "#4e5c51");
        yTick.setAttribute("font-size", "28");
        yTick.textContent = gy.toFixed(2);
        bifurcationSvg.appendChild(yTick);

        var yTickMark = makeSvgEl("line");
        yTickMark.setAttribute("x1", margin.left - 7);
        yTickMark.setAttribute("x2", margin.left);
        yTickMark.setAttribute("y1", yScale(gy));
        yTickMark.setAttribute("y2", yScale(gy));
        yTickMark.setAttribute("stroke", "#425044");
        yTickMark.setAttribute("stroke-width", "1.5");
        bifurcationSvg.appendChild(yTickMark);
      }

      var axisX = makeSvgEl("line");
      axisX.setAttribute("x1", margin.left);
      axisX.setAttribute("x2", margin.left + innerWidth);
      axisX.setAttribute("y1", margin.top + innerHeight);
      axisX.setAttribute("y2", margin.top + innerHeight);
      axisX.setAttribute("stroke", "#425044");
      axisX.setAttribute("stroke-width", "2");
      bifurcationSvg.appendChild(axisX);

      var axisY = makeSvgEl("line");
      axisY.setAttribute("x1", margin.left);
      axisY.setAttribute("x2", margin.left);
      axisY.setAttribute("y1", margin.top);
      axisY.setAttribute("y2", margin.top + innerHeight);
      axisY.setAttribute("stroke", "#425044");
      axisY.setAttribute("stroke-width", "2");
      bifurcationSvg.appendChild(axisY);

      var xLabel = makeSvgEl("text");
      xLabel.setAttribute("x", margin.left + innerWidth / 2);
      xLabel.setAttribute("y", height - 6);
      xLabel.setAttribute("text-anchor", "middle");
      xLabel.setAttribute("fill", "#1f2c22");
      xLabel.setAttribute("font-size", "28");
      xLabel.textContent = "slope γ [°]";
      bifurcationSvg.appendChild(xLabel);

      var yLabel = makeSvgEl("text");
      yLabel.setAttribute("x", "20");
      yLabel.setAttribute("y", margin.top + innerHeight / 2);
      yLabel.setAttribute("transform", "rotate(-90 20 " + (margin.top + innerHeight / 2) + ")");
      yLabel.setAttribute("text-anchor", "middle");
      yLabel.setAttribute("fill", "#1f2c22");
      yLabel.setAttribute("font-size", "28");
      yLabel.textContent = "period T [t₀]";
      bifurcationSvg.appendChild(yLabel);

      [
        { name: "short", color: "#2f6f3e", data: state.continuation ? state.continuation.short : [] },
        { name: "long", color: "#7b4aa1", data: state.continuation ? state.continuation.long : [] }
      ].forEach(function (branch) {
        if (!branch.data || !branch.data.length) return;
        [-1, 1].forEach(function (mirrorSign) {
          var d = "";
          branch.data.forEach(function (point) {
            var gammaDeg = clamp(mirrorSign * point.gammaDeg, -gammaMax, gammaMax);
            var pointX = xScale(gammaDeg);
            var pointY = yScale(point.T);
            d += (d ? " L " : "M ") + pointX.toFixed(2) + " " + pointY.toFixed(2);
            state.bifurcationHitPoints.push({
              x: pointX,
              y: pointY,
              gammaDeg: mirrorSign * point.gammaDeg,
              branch: branch.name,
              T: point.T
            });
          });
          var path = makeSvgEl("path");
          path.setAttribute("d", d);
          path.setAttribute("fill", "none");
          path.setAttribute("stroke", branch.color);
          path.setAttribute("stroke-width", state.branch === branch.name ? "6" : "5");
          path.setAttribute("stroke-linecap", "round");
          path.setAttribute("stroke-linejoin", "round");
          path.setAttribute("opacity", state.branch === branch.name ? "1" : "0.75");
          bifurcationSvg.appendChild(path);
        });

        var markerPoint = closestBranchPoint(branch.data, Math.abs(state.gammaDeg)) || branch.data[branch.data.length - 1];
        var markerSign = state.gammaDeg < 0 ? -1 : 1;
        var branchStartPoint = branch.data[0] || markerPoint;
        var markerGamma = state.branch === branch.name && state.currentOrbit
          ? (state.currentOrbit.isStandstill ? 0 : state.currentOrbit.gammaDeg)
          : markerSign * markerPoint.gammaDeg;
        var markerT = state.branch === branch.name && state.currentOrbit
          ? (state.currentOrbit.isStandstill ? branchStartPoint.T : state.currentOrbit.T)
          : markerPoint.T;
        var marker = makeSvgEl("circle");
        marker.setAttribute("cx", xScale(clamp(markerGamma, -gammaMax, gammaMax)));
        marker.setAttribute("cy", yScale(markerT));
        marker.setAttribute("r", "11");
        marker.setAttribute("fill", state.branch === branch.name ? branch.color : "#ffffff");
        marker.setAttribute("stroke", branch.color);
        marker.setAttribute("stroke-width", "2.5");
        bifurcationSvg.appendChild(marker);
      });

      updateBifurcationCursor(state.lastBifurcationPointer);
    }

    function updateOrbitOnly() {
      renderControlValues();
      rebuildContinuationIfNeeded();
      state.currentOrbit = solveOrbitForGamma(state.gammaDeg, state.branch, currentParams(), state.continuation);
      state.groundMarkerLastAnchorWorldX = null;
      drawBifurcationDiagram();
      drawWalker();
    }

    function scheduleOrbitUpdate() {
      if (state.orbitUpdateScheduled) return;
      state.orbitUpdateScheduled = true;
      requestAnimationFrame(function () {
        state.orbitUpdateScheduled = false;
        updateOrbitOnly();
      });
    }

    function updateParametersWithLoading() {
      var token;
      syncStateFromControls();
      renderControlValues();
      setLoading(true, "compute new bifurcation diagram ...");
      token = ++state.loadingToken;
      requestAnimationFrame(function () {
        requestAnimationFrame(function () {
          if (token !== state.loadingToken) return;
          rebuildContinuationIfNeeded();
          state.currentOrbit = solveOrbitForGamma(state.gammaDeg, state.branch, currentParams(), state.continuation);
          state.groundMarkerLastAnchorWorldX = null;
          drawBifurcationDiagram();
          drawWalker();
          setLoading(false);
        });
      });
    }

    function canvasPointFromEvent(event) {
      state.lastCanvasPointer = eventToElementCoordinates(event, canvas, canvas.width, canvas.height);
      return state.lastCanvasPointer;
    }

    function bifurcationPointFromEvent(event) {
      state.lastBifurcationPointer = eventToElementCoordinates(event, bifurcationSvg, bifurcationSize.width, bifurcationSize.height);
      return state.lastBifurcationPointer;
    }

    function getSlopeDragTarget(point) {
      var interaction = state.slopeInteraction;
      var leftDistance;
      var rightDistance;
      if (!interaction || !point) return null;
      leftDistance = pointDistance(point, interaction.leftHandle);
      rightDistance = pointDistance(point, interaction.rightHandle);
      if (leftDistance <= interaction.handleRadius + 8 && leftDistance <= rightDistance) return "left";
      if (rightDistance <= interaction.handleRadius + 8) return "right";
      if (distanceToSegment(point, interaction.leftHandle, interaction.rightHandle) <= interaction.lineRadius) {
        return point.x <= interaction.midpoint.x ? "left" : "right";
      }
      return null;
    }

    function updateCanvasCursor(point) {
      if (state.slopeDrag) {
        canvas.style.cursor = "grabbing";
        return;
      }
      canvas.style.cursor = getSlopeDragTarget(point) ? "pointer" : "default";
    }

    function updateGammaFromSlopeDrag(point) {
      var dy;
      var gammaRad;
      var nextGammaDeg;
      if (!state.slopeDrag) return;
      dy = point.y - state.slopeDrag.midpoint.y;
      gammaRad = state.slopeDrag.side === "right"
        ? Math.atan2(dy, state.slopeDrag.halfWidth)
        : Math.atan2(-dy, state.slopeDrag.halfWidth);
      nextGammaDeg = clamp(gammaRad * 180 / Math.PI, gammaRange.min, gammaRange.max);
      nextGammaDeg = normalizeGammaDeg(nextGammaDeg, gammaRange.min, gammaRange.max);
      if (Math.abs(nextGammaDeg - state.gammaDeg) < 1e-4) return;
      state.gammaDeg = nextGammaDeg;
      renderControlValues();
      scheduleOrbitUpdate();
    }

    function findNearestBifurcationHit(point) {
      var best = null;
      var bestDistance = Infinity;
      state.bifurcationHitPoints.forEach(function (candidate) {
        var distance = pointDistance(point, candidate);
        if (distance < bestDistance) {
          best = candidate;
          bestDistance = distance;
        }
      });
      return bestDistance <= 18 ? best : null;
    }

    function updateBifurcationCursor(point) {
      bifurcationSvg.style.cursor = point && findNearestBifurcationHit(point) ? "pointer" : "default";
    }

    controls.a.addEventListener("input", updateParametersWithLoading);
    controls.m.addEventListener("input", updateParametersWithLoading);
    controls.branch.forEach(function (input) {
      input.addEventListener("change", function () {
        syncStateFromControls();
        updateOrbitOnly();
      });
    });

    canvas.addEventListener("pointerdown", function (event) {
      var point;
      var targetSide;
      if (interactive.classList.contains("is-loading")) return;
      point = canvasPointFromEvent(event);
      targetSide = getSlopeDragTarget(point);
      if (!targetSide) {
        updateCanvasCursor(point);
        return;
      }
      event.preventDefault();
      state.slopeDrag = {
        pointerId: event.pointerId,
        side: targetSide,
        midpoint: {
          x: state.slopeInteraction.midpoint.x,
          y: state.slopeInteraction.midpoint.y
        },
        halfWidth: Math.max(Math.abs(state.slopeInteraction.rightHandle.x - state.slopeInteraction.leftHandle.x) * 0.5, 1)
      };
      if (canvas.setPointerCapture) canvas.setPointerCapture(event.pointerId);
      updateCanvasCursor(point);
      updateGammaFromSlopeDrag(point);
    });

    canvas.addEventListener("pointermove", function (event) {
      var point = canvasPointFromEvent(event);
      if (state.slopeDrag && state.slopeDrag.pointerId === event.pointerId) {
        updateGammaFromSlopeDrag(point);
        return;
      }
      updateCanvasCursor(point);
    });

    function endSlopeDrag(event) {
      if (state.slopeDrag && state.slopeDrag.pointerId === event.pointerId) {
        if (canvas.releasePointerCapture) canvas.releasePointerCapture(event.pointerId);
        state.slopeDrag = null;
      }
      if (typeof event.clientX === "number" && typeof event.clientY === "number") {
        updateCanvasCursor(canvasPointFromEvent(event));
      } else {
        canvas.style.cursor = "default";
      }
    }

    canvas.addEventListener("pointerup", endSlopeDrag);
    canvas.addEventListener("pointercancel", endSlopeDrag);
    canvas.addEventListener("pointerleave", function () {
      if (state.slopeDrag) return;
      state.lastCanvasPointer = null;
      canvas.style.cursor = "default";
    });

    bifurcationSvg.addEventListener("pointermove", function (event) {
      updateBifurcationCursor(bifurcationPointFromEvent(event));
    });

    bifurcationSvg.addEventListener("pointerleave", function () {
      state.lastBifurcationPointer = null;
      bifurcationSvg.style.cursor = "default";
    });

    bifurcationSvg.addEventListener("click", function (event) {
      var hit = findNearestBifurcationHit(bifurcationPointFromEvent(event));
      if (!hit) return;
      state.branch = hit.branch;
      state.gammaDeg = normalizeGammaDeg(hit.gammaDeg, gammaRange.min, gammaRange.max);
      updateOrbitOnly();
    });

    function animate(timestamp) {
      if (state.lastFrame == null) state.lastFrame = timestamp;
      var dt = Math.min(0.05, (timestamp - state.lastFrame) / 1000);
      state.lastFrame = timestamp;
      state.time += dt;
      if (state.currentOrbit) drawWalker();
      requestAnimationFrame(animate);
    }

    syncStateFromControls();
    updateOrbitOnly();
    requestAnimationFrame(animate);
  });
}());

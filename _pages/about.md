---
permalink: /
title: ""
author_profile: true
redirect_from: 
  - /about/
  - /about.html
---

Hi, I’m Max. I grew up just south of Stuttgart, Germany. From early on, I enjoyed tackling intricate projects—from carpentry around the house to transcribing and arranging orchestral music. What fascinated me was how combining many pieces could create something that feels simple, yet beautiful and innovative.

Today, I apply mathematics to real-world systems such as robotics, bringing together ideas from different areas to develop elegant and practical solutions. (<a href="/research/">Watch my research video to learn more.</a>)

Outside the lab, I’m a passionate trombone player. I regularly perform with Stuttgart’s salsa ensemble Caballo Negro <a href="https://www.youtube.com/@caballo-negro" target="_blank" rel="noopener" aria-label="Caballo Negro on YouTube"><i class="fab fa-youtube" aria-hidden="true"></i></a> and also lead my own jazz band, Youngblood Dixies <a href="https://www.youtube.com/@youngblooddixies" target="_blank" rel="noopener" aria-label="Youngblood Dixies on YouTube"><i class="fab fa-youtube" aria-hidden="true"></i></a>.

Professional Journey and PhD Research
======
I am currently a postdoctoral researcher at the University of Stuttgart, working on the **co-design of hardware and behavior for legged robotic systems** in the newly founded robotics group at the <a href="https://www.iams.uni-stuttgart.de/" target="_blank" rel="noopener">Institute for Adaptive Mechanical Systems (IAMS)</a>. The group is located within Cyber Valley at the Max Planck Institute for Intelligent Systems.

Previously, during my PhD, I worked on new mathematical concepts for energy-efficient locomotion. Below is a short video introducing my PhD research.
<video id="thesis-video" width="640" height="360" controls playsinline preload="metadata" style="border: 1px solid #b3b3b3;">
  <source src="/files/Raff_Maximilian_ThesisVideo_web.mp4" type="video/mp4">
  <track id="thesis-video-captions" kind="subtitles" srclang="en" label="English" src="/files/Raff_Maximilian_ThesisVideo.vtt" default>
  Your browser does not support the video tag.
</video>
<p id="thesis-video-fallback" style="display:none;">
  If embedded playback fails, open the video directly:
  <a href="/files/Raff_Maximilian_ThesisVideo_web.mp4" target="_blank" rel="noopener">Raff_Maximilian_ThesisVideo_web.mp4</a>
</p>

<script>
(function () {
  var video = document.getElementById("thesis-video");
  if (!video) return;

  function updateCaptionMode() {
    var mode = (video.muted || video.volume === 0) ? "showing" : "hidden";
    if (!video.textTracks) return;
    for (var i = 0; i < video.textTracks.length; i++) {
      var track = video.textTracks[i];
      if (track.kind === "subtitles" || track.kind === "captions") {
        try {
          track.mode = mode;
        } catch (e) {
          // Ignore per-browser track mode errors.
        }
      }
    }
  }

  video.addEventListener("loadedmetadata", updateCaptionMode);
  video.addEventListener("loadeddata", updateCaptionMode);
  video.addEventListener("volumechange", updateCaptionMode);
  video.addEventListener("error", function () {
    var fallback = document.getElementById("thesis-video-fallback");
    if (fallback) fallback.style.display = "block";
  });
  updateCaptionMode();
})();
</script>

My PhD, mentored by <a href="https://scholar.google.de/citations?user=tx_oGsMAAAAJ&hl=de&oi=sra" target="_blank" rel="noopener">C. David Remy</a>, resulted in the thesis <a href="https://elib.uni-stuttgart.de/items/3c4b9711-0be2-41a5-ab28-8f86a3d06b49" target="_blank" rel="noopener">“Continuation and Optimization of Gaits and Other Non-Smooth Orbits”</a>, awarded the highest distinction (summa cum laude).

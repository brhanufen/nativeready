/* NativeReady cosmetic motion only.
 *
 * This file is ADDITIVE and isolated. It never calls /predict or /feedback, never
 * reads or writes the sequence, and never changes the value script.js sends to the
 * server (feedback posts currentPrediction.suitability_score, an internal JS value,
 * not the DOM text this file animates). It only: (1) reveals .reveal elements on
 * scroll, and (2) animates the *displayed* score number. Fully gated by
 * prefers-reduced-motion; if reduced motion is requested, this file does nothing.
 */
(function () {
  "use strict";

  var reduce = window.matchMedia &&
    window.matchMedia("(prefers-reduced-motion: reduce)").matches;
  if (reduce) return; // respect the user's preference: add no cosmetic motion

  document.documentElement.classList.add("motion");

  // --- P2: staggered scroll-reveal via IntersectionObserver (works in all browsers) ---
  function initReveals() {
    var els = Array.prototype.slice.call(document.querySelectorAll(".reveal"));
    if (!els.length) return;
    if (!("IntersectionObserver" in window)) {
      els.forEach(function (el) { el.classList.add("is-visible"); });
      return;
    }
    var io = new IntersectionObserver(function (entries) {
      entries.forEach(function (entry) {
        if (!entry.isIntersecting) return;
        var el = entry.target;
        var group = Array.prototype.slice.call(
          el.parentNode.querySelectorAll(".reveal"));
        var idx = group.indexOf(el);
        el.style.transitionDelay = Math.min((idx < 0 ? 0 : idx) * 80, 320) + "ms";
        el.classList.add("is-visible");
        io.unobserve(el);
      });
    }, { threshold: 0.12, rootMargin: "0px 0px -8% 0px" });
    els.forEach(function (el) { io.observe(el); });
  }

  // --- P3: score count-up (display-only) ---
  function initCountUp() {
    var scoreEl = document.getElementById("result-score");
    if (!scoreEl || !("MutationObserver" in window) ||
        typeof window.requestAnimationFrame !== "function") return;

    var shown = null, rafId = null;
    var cfg = { childList: true, characterData: true, subtree: true };

    var mo = new MutationObserver(function () {
      var target = parseInt((scoreEl.textContent || "").trim(), 10);
      if (isNaN(target) || target === shown) return;
      var from = (shown === null) ? 0 : shown;
      shown = target;
      if (rafId) cancelAnimationFrame(rafId);
      mo.disconnect(); // do not observe our own intermediate writes
      var start = null, dur = 700;
      function step(ts) {
        if (start === null) start = ts;
        var p = Math.min(1, (ts - start) / dur);
        var eased = 1 - Math.pow(1 - p, 3); // easeOutCubic
        scoreEl.textContent = String(Math.round(from + (target - from) * eased));
        if (p < 1) {
          rafId = requestAnimationFrame(step);
        } else {
          scoreEl.textContent = String(target); // land exactly on the real value
          mo.observe(scoreEl, cfg);
        }
      }
      rafId = requestAnimationFrame(step);
    });
    mo.observe(scoreEl, cfg);
  }

  // --- On-brand hero animation: a subtle drifting protein-sequence field ---
  function initSeqField() {
    var canvas = document.querySelector(".seq-field");
    if (!canvas || !canvas.getContext || typeof window.requestAnimationFrame !== "function") return;
    var ctx = canvas.getContext("2d");
    var AAs = "ACDEFGHIKLMNPQRSTVWY"; // the 20 standard amino acids the tool ingests
    var dpr = Math.min(window.devicePixelRatio || 1, 2);
    var W = 0, H = 0, parts = [], raf = null, last = 0;

    function resize() {
      var r = canvas.getBoundingClientRect();
      W = Math.max(1, r.width); H = Math.max(1, r.height);
      canvas.width = Math.floor(W * dpr); canvas.height = Math.floor(H * dpr);
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    }
    function make(seed) {
      return {
        x: Math.random() * W,
        y: seed ? Math.random() * H : H + 24,
        ch: AAs.charAt((Math.random() * AAs.length) | 0),
        size: 11 + Math.random() * 8,
        vy: 7 + Math.random() * 12,          // px/sec, slow rise
        vx: (Math.random() - 0.5) * 5,
        life: 0, ttl: 6 + Math.random() * 8
      };
    }
    function seed() {
      var n = Math.max(26, Math.min(64, Math.floor(W / 26)));
      parts = []; for (var i = 0; i < n; i++) parts.push(make(true));
    }
    function frame(ts) {
      if (!last) last = ts;
      var dt = Math.min(0.05, (ts - last) / 1000); last = ts;
      ctx.clearRect(0, 0, W, H);
      ctx.textAlign = "center"; ctx.textBaseline = "middle";
      for (var i = 0; i < parts.length; i++) {
        var p = parts[i];
        p.life += dt; p.y -= p.vy * dt; p.x += p.vx * dt;
        var t = p.life / p.ttl;
        var fade = t < 0.2 ? t / 0.2 : (t > 0.8 ? (1 - t) / 0.2 : 1);
        var a = Math.max(0, fade) * 0.30;
        ctx.font = p.size.toFixed(1) + "px ui-monospace, SFMono-Regular, Menlo, monospace";
        ctx.fillStyle = "rgba(63, 224, 210, " + a.toFixed(3) + ")";
        ctx.fillText(p.ch, p.x, p.y);
        if (p.life >= p.ttl || p.y < -24) parts[i] = make(false);
      }
      raf = requestAnimationFrame(frame);
    }
    function start() { if (!raf) { last = 0; raf = requestAnimationFrame(frame); } }
    function stop() { if (raf) { cancelAnimationFrame(raf); raf = null; } }

    resize(); seed();
    var rt;
    window.addEventListener("resize", function () {
      clearTimeout(rt); rt = setTimeout(function () { resize(); seed(); }, 150);
    });
    document.addEventListener("visibilitychange", function () {
      if (document.hidden) stop(); else start();
    });
    start();
  }

  // --- Hero video: play only when motion is allowed (this file already returned
  //     early under reduced motion, so reduced-motion users keep the static poster). ---
  function initHeroVideo() {
    var v = document.querySelector(".hero-video");
    if (!v) return false;
    try {
      v.muted = true;
      var pr = v.play();
      if (pr && typeof pr.catch === "function") pr.catch(function () {});
    } catch (e) {}
    return true;
  }

  // --- Dynamic rotating note under the input (honest one-liners, fade cycle) ---
  function initRotatingNote() {
    var el = document.querySelector(".rn-text");
    if (!el) return;
    var notes = [
      "Membrane proteins are flagged for nanodiscs / amphipols.",
      "Glycoproteins are flagged for PNGase F deglycosylation.",
      "Scores are a first-pass prioritizer, not a verdict.",
      "Open 635-protein benchmark · ESM-2 + BioPython features.",
      "Risk panel: molecular weight, pI, membrane topology, glycosylation."
    ];
    var i = 0;
    setInterval(function () {
      el.style.opacity = "0";
      setTimeout(function () {
        i = (i + 1) % notes.length;
        el.textContent = notes[i];
        el.style.opacity = "1";
      }, 420);
    }, 4200);
  }

  // On brand hero motion: a faint, living native mass spectrum behind the headline.
  // Peaks breathe, a slow highlight sweeps left to right. Color comes from the theme
  // accent so it adapts to light and dark. Purely decorative canvas, no DOM contract.
  function initHeroSpectrum() {
    var canvas = document.querySelector(".hero-spectrum");
    if (!canvas || !canvas.getContext || typeof window.requestAnimationFrame !== "function") return;
    var ctx = canvas.getContext("2d");
    var dpr = Math.min(window.devicePixelRatio || 1, 2);
    var W = 0, H = 0, raf = null, t0 = 0, peaks = [], N = 0;

    function accent() {
      var c = getComputedStyle(document.documentElement).getPropertyValue("--accent").trim();
      return c || "#10935F";
    }
    function build() {
      N = Math.max(30, Math.min(72, Math.floor(W / 20)));
      peaks = [];
      for (var i = 0; i < N; i++) {
        var x = i / (N - 1);
        // two gaussian humps read as the charge state / glycoform envelope of a native spectrum
        var env = Math.exp(-Math.pow((x - 0.40) / 0.17, 2)) * 0.95
                + Math.exp(-Math.pow((x - 0.66) / 0.11, 2)) * 0.60 + 0.05;
        peaks.push({ x: x, base: env * (0.34 + Math.random() * 0.66),
                     ph: Math.random() * 6.283, sp: 0.45 + Math.random() * 0.9 });
      }
    }
    function resize() {
      var r = canvas.getBoundingClientRect();
      W = Math.max(1, r.width); H = Math.max(1, r.height);
      canvas.width = Math.floor(W * dpr); canvas.height = Math.floor(H * dpr);
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0); build();
    }
    function frame(ts) {
      if (!t0) t0 = ts;
      var t = (ts - t0) / 1000, col = accent();
      var baseY = H * 0.98, maxH = H * 0.80;
      var scan = (t * 0.11) % 1.25 - 0.12; // slow highlight sweep across the envelope
      ctx.clearRect(0, 0, W, H);
      for (var i = 0; i < N; i++) {
        var p = peaks[i];
        var breath = 0.80 + 0.20 * Math.sin(t * p.sp + p.ph);
        var h = p.base * maxH * breath, x = p.x * W;
        var glow = Math.max(0, 1 - Math.min(1, Math.abs(p.x - scan) / 0.10));
        ctx.strokeStyle = col;
        ctx.globalAlpha = 0.13 + p.base * 0.10 + glow * 0.22;
        ctx.lineWidth = Math.max(1.1, (W / N) * 0.13);
        ctx.beginPath(); ctx.moveTo(x, baseY); ctx.lineTo(x, baseY - h); ctx.stroke();
        ctx.globalAlpha = 0.22 * breath + glow * 0.40; ctx.fillStyle = col;
        ctx.beginPath(); ctx.arc(x, baseY - h, 1.5 + glow * 1.2, 0, 6.283); ctx.fill();
      }
      ctx.globalAlpha = 0.14; ctx.strokeStyle = col; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(0, baseY); ctx.lineTo(W, baseY); ctx.stroke();
      ctx.globalAlpha = 1;
      raf = requestAnimationFrame(frame);
    }
    function start() { if (!raf) { t0 = 0; raf = requestAnimationFrame(frame); } }
    function stop() { if (raf) { cancelAnimationFrame(raf); raf = null; } }
    resize();
    var rt;
    window.addEventListener("resize", function () { clearTimeout(rt); rt = setTimeout(resize, 150); });
    document.addEventListener("visibilitychange", function () { if (document.hidden) stop(); else start(); });
    start();
  }

  // Auto-playing product demo: type a sequence, predict, reveal the score, loop.
  // Purely cosmetic; the real predictor is the live tool above. Under reduced
  // motion this file returns early, so the demo stays in its static filled state.
  function initDemo() {
    var root = document.querySelector(".demo");
    if (!root) return;
    var input  = root.querySelector(".demo-input"),
        typed  = root.querySelector(".demo-typed"),
        btn    = root.querySelector(".demo-btn"),
        result = root.querySelector(".demo-result"),
        num    = root.querySelector(".demo-num"),
        fill   = root.querySelector(".demo-bar-fill"),
        pills  = Array.prototype.slice.call(root.querySelectorAll(".demo-pill"));
    if (!typed || !btn || !result || !num || !fill) return;
    var SEQ = "MQIFVKTLTGKTITLEVEPSDTIENVK...";
    var timers = [];
    function later(fn, ms) { timers.push(setTimeout(fn, ms)); }
    function clearAll() { timers.forEach(clearTimeout); timers = []; }
    function reset() {
      clearAll();
      typed.textContent = "";
      if (input) input.classList.add("active");
      result.hidden = true;
      num.textContent = "0";
      fill.style.width = "0%";
      pills.forEach(function (p) { p.style.opacity = "0"; });
      btn.classList.remove("loading");
    }
    function type(i) {
      if (i > SEQ.length) { afterType(); return; }
      typed.textContent = SEQ.slice(0, i);
      later(function () { type(i + 1); }, 24 + Math.random() * 34);
    }
    function countUp(to) {
      var start = null, dur = 800;
      function step(ts) {
        if (start === null) start = ts;
        var p = Math.min(1, (ts - start) / dur);
        num.textContent = String(Math.round(to * (1 - Math.pow(1 - p, 3))));
        if (p < 1) requestAnimationFrame(step); else num.textContent = String(to);
      }
      requestAnimationFrame(step);
    }
    function afterType() {
      if (input) input.classList.remove("active");
      later(function () {
        btn.classList.add("loading");
        later(function () {
          btn.classList.remove("loading");
          result.hidden = false;
          requestAnimationFrame(function () { fill.style.width = "96%"; });
          countUp(96);
          pills.forEach(function (p, idx) { later(function () { p.style.opacity = "1"; }, 180 + idx * 150); });
          later(loop, 5000);
        }, 1100);
      }, 600);
    }
    function loop() { reset(); later(function () { type(1); }, 650); }
    document.addEventListener("visibilitychange", function () {
      if (document.hidden) clearAll(); else loop();
    });
    loop();
  }

  function init() {
    var hasVideo = initHeroVideo();
    initReveals();
    initCountUp();
    initRotatingNote();
    initHeroSpectrum();
    initDemo();
    if (!hasVideo) initSeqField(); // canvas field is the fallback when no video
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();

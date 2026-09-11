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

  function init() {
    var hasVideo = initHeroVideo();
    initReveals();
    initCountUp();
    initRotatingNote();
    if (!hasVideo) initSeqField(); // canvas field is the fallback when no video
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();

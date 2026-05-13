// NativeReady frontend
// Plain JS, no build step. Configure API_BASE for your deployment.

const API_BASE = window.location.hostname === 'localhost'
  ? 'http://localhost:8000'
  : 'https://nativeready-production.up.railway.app';

// Real UniProt sequences for sample buttons. All canonical, verifiable.
const SAMPLES = {
  carbonic_anhydrase: {
    label: 'Carbonic anhydrase 2 (P00918)',
    seq: '>sp|P00918|CAH2_HUMAN Carbonic anhydrase 2 OS=Homo sapiens\n' +
         'MSHHWGYGKHNGPEHWHKDFPIAKGERQSPVDIDTHTAKYDPSLKPLSVSYDQATSLRILN\n' +
         'NGHAFNVEFDDSQDKAVLKGGPLDGTYRLIQFHFHWGSLDGQGSEHTVDKKKYAAELHLVH\n' +
         'WNTKYGDFGKAVQQPDGLAVLGIFLKVGSAKPGLQKVVDVLDSIKTKGKSADFTNFDPRGL\n' +
         'LPESLDYWTYPGSLTTPPLLECVTWIVLKEPISVSSEQVLKFRKLNFNGEGEPEELMVDNW\n' +
         'RPAQPLKNRQIKASFK',
  },
  ubiquitin: {
    label: 'Ubiquitin (P0CG48)',
    seq: '>sp|P0CG48|UBC_HUMAN Polyubiquitin-C OS=Homo sapiens\n' +
         'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNI\n' +
         'QKESTLHLVLRLRGG',
  },
  lysozyme: {
    label: 'Lysozyme C, chicken (P00698)',
    seq: '>sp|P00698|LYSC_CHICK Lysozyme C OS=Gallus gallus\n' +
         'MRSLLILVLCFLPLAALGKVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQAT\n' +
         'NRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNG\n' +
         'MNAWVAWRNRCKGTDVQAWIRGCRL',
  },
  insulin: {
    label: 'Insulin chain B (small peptide)',
    seq: '>insulin_B_chain Human insulin B chain (30 aa)\n' +
         'FVNQHLCGSHLVEALYLVCGERGFFYTPKT',
  },
  mucin: {
    label: 'Mucin-1 (Q9Y6Y8 fragment, hard case)',
    seq: '>mucin_repeat Tandem repeat region of human MUC1\n' +
         'GSTAPPAHGVTSAPDTRPAPGSTAPPAHGVTSAPDTRPAPGSTAPPAHGVTSAPDTRPAP' +
         'GSTAPPAHGVTSAPDTRPAPGSTAPPAHGVTSAPDTRPAPGSTAPPAHGVTSAPDTRPAP',
  },
};

// ---- DOM refs
const form        = document.getElementById('predict-form');
const textarea    = document.getElementById('sequence');
const submitBtn   = document.getElementById('submit-btn');
const sampleBtns  = document.querySelectorAll('.sample-btn');
const feedbackCard   = document.getElementById('feedback-card');
const feedbackThanks = document.getElementById('feedback-thanks');
const feedbackBtns   = document.querySelectorAll('.feedback-btn');

// Track current prediction for feedback context
let currentPrediction = null;
let currentSequence = null;
const resetBtn    = document.getElementById('reset-btn');
const errorBox    = document.getElementById('error-box');
const resultsBox  = document.getElementById('results');
const scoreEl     = document.getElementById('result-score');
const labelEl     = document.getElementById('result-label');
const ciEl        = document.getElementById('result-ci');
const barEl       = document.getElementById('result-bar');
const risksEl     = document.getElementById('result-risks');
const recsEl      = document.getElementById('result-recs');
const modelVerEl  = document.getElementById('model-version');

// ---- Utilities

// Strip FASTA headers and whitespace, uppercase, keep only valid residue chars.
function cleanSequence(raw) {
  if (!raw) return '';
  const lines = raw.split(/\r?\n/);
  const body = lines.filter(l => !l.startsWith('>')).join('');
  return body.replace(/\s+/g, '').toUpperCase();
}

function isPlausibleProtein(seq) {
  if (seq.length < 20) return false;
  // Standard amino acids plus ambiguity codes
  return /^[ACDEFGHIKLMNPQRSTVWYBXZJUO*]+$/.test(seq);
}

function scoreClass(score) {
  if (score >= 70) return 'good';
  if (score >= 40) return 'mid';
  return 'bad';
}

function showError(msg) {
  errorBox.innerHTML = '<strong>Something went wrong.</strong>' + escapeHtml(msg);
  errorBox.hidden = false;
  resultsBox.hidden = true;
}

function clearError() {
  errorBox.hidden = true;
  errorBox.textContent = '';
}

function escapeHtml(s) {
  return String(s)
    .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;').replace(/'/g, '&#39;');
}

function setLoading(isLoading) {
  if (isLoading) {
    submitBtn.disabled = true;
    submitBtn.classList.add('is-loading');
    submitBtn.querySelector('.btn-label').textContent = 'Predicting...';
  } else {
    submitBtn.disabled = false;
    submitBtn.classList.remove('is-loading');
    submitBtn.querySelector('.btn-label').textContent = 'Predict suitability';
  }
}

// ---- Render

function renderResults(data) {
  const score = clamp(Number(data.suitability_score) || 0, 0, 100);
  const cls = scoreClass(score);

  scoreEl.textContent = Math.round(score);
  scoreEl.className = 'score-num ' + cls;

  labelEl.textContent = data.suitability_label || labelFromScore(score);

  // OOD warning if present
  let oodEl = document.getElementById('ood-warning');
  if (oodEl) oodEl.remove();
  if (data.warning === 'out_of_distribution') {
    oodEl = document.createElement('div');
    oodEl.id = 'ood-warning';
    oodEl.className = 'ood-warning';
    oodEl.innerHTML = '<strong>Low confidence — out of distribution.</strong> ' +
      escapeHtml(data.warning_message || 'This sequence is unusual relative to training data. Treat the score with extra caution.');
    resultsBox.insertBefore(oodEl, resultsBox.firstChild);
  }

  if (data.confidence_interval &&
      data.confidence_interval.lower != null &&
      data.confidence_interval.upper != null) {
    const lo = Math.round(data.confidence_interval.lower);
    const hi = Math.round(data.confidence_interval.upper);
    ciEl.textContent = '95% confidence interval: ' + lo + '–' + hi;
    ciEl.hidden = false;
  } else {
    ciEl.textContent = '';
    ciEl.hidden = true;
  }

  // animate bar
  barEl.className = 'bar-fill ' + cls;
  // force reflow so transition runs from 0
  barEl.style.width = '0%';
  // eslint-disable-next-line no-unused-expressions
  barEl.offsetWidth;
  barEl.style.width = score + '%';

  // risks
  risksEl.innerHTML = '';
  const risks = Array.isArray(data.risk_factors) ? data.risk_factors : [];
  if (risks.length === 0) {
    risksEl.innerHTML = '<div class="risk-card"><div class="rname">No risk factors reported.</div></div>';
  } else {
    risks.forEach(r => {
      const lvl = String(r.risk_level || 'unknown').toLowerCase();
      const card = document.createElement('div');
      card.className = 'risk-card';
      card.innerHTML =
        '<div class="rname">' + escapeHtml(r.name || 'Factor') + '</div>' +
        '<div class="rval">' + escapeHtml(r.value != null ? r.value : '—') + '</div>' +
        '<span class="risk-pill ' + escapeHtml(lvl) + '">' + escapeHtml(lvl) + ' risk</span>';
      risksEl.appendChild(card);
    });
  }

  // recommendations
  recsEl.innerHTML = '';
  const recs = Array.isArray(data.recommendations) ? data.recommendations : [];
  if (recs.length === 0) {
    const li = document.createElement('li');
    li.textContent = 'No specific recommendations returned.';
    recsEl.appendChild(li);
  } else {
    recs.forEach(rec => {
      const li = document.createElement('li');
      li.textContent = rec;
      recsEl.appendChild(li);
    });
  }

  modelVerEl.textContent = data.model_version
    ? 'Model: ' + data.model_version
    : '';

  resultsBox.hidden = false;
  // smooth scroll into view
  setTimeout(() => {
    resultsBox.scrollIntoView({ behavior: 'smooth', block: 'start' });
  }, 60);
}

function clamp(n, lo, hi) { return Math.max(lo, Math.min(hi, n)); }

function labelFromScore(s) {
  if (s >= 80) return 'Excellent';
  if (s >= 65) return 'Good';
  if (s >= 45) return 'Marginal';
  if (s >= 25) return 'Poor';
  return 'Very poor';
}

// ---- API

async function predict(sequence) {
  let res;
  try {
    res = await fetch(API_BASE + '/predict', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ sequence: sequence }),
    });
  } catch (e) {
    throw new Error('Could not reach the prediction server. Check your connection and try again.');
  }

  let payload = null;
  try { payload = await res.json(); } catch (_) { /* ignore */ }

  if (!res.ok) {
    const serverMsg =
      (payload && (payload.detail || payload.error || payload.message)) ||
      ('Server returned status ' + res.status + '.');
    throw new Error(serverMsg);
  }

  if (!payload || typeof payload !== 'object') {
    throw new Error('Server returned an empty or invalid response.');
  }

  return payload;
}

// ---- Events

form.addEventListener('submit', async (ev) => {
  ev.preventDefault();
  clearError();

  const cleaned = cleanSequence(textarea.value);
  if (!cleaned) {
    showError('Please paste a protein sequence first.');
    return;
  }
  if (!isPlausibleProtein(cleaned)) {
    showError(
      'That doesn\'t look like a protein sequence. Use single-letter amino acid codes (at least 20 residues).'
    );
    return;
  }

  setLoading(true);
  resultsBox.hidden = true;

  try {
    const data = await predict(cleaned);
    currentPrediction = data;
    currentSequence = cleaned;
    // Reset feedback UI for the new prediction
    feedbackCard.hidden = false;
    feedbackThanks.hidden = true;
    feedbackBtns.forEach(b => b.disabled = false);
    renderResults(data);
  } catch (err) {
    showError(err.message || 'Unknown error.');
  } finally {
    setLoading(false);
  }
});

// ---- Feedback handler
feedbackBtns.forEach(btn => {
  btn.addEventListener('click', async () => {
    if (!currentPrediction || !currentSequence) return;
    const outcome = btn.dataset.outcome;
    feedbackBtns.forEach(b => b.disabled = true);
    try {
      const emailEl = document.getElementById('feedback-email');
      const noteEl  = document.getElementById('feedback-note');
      const email = emailEl ? emailEl.value.trim() : '';
      const note  = noteEl ? noteEl.value.trim() : '';
      const resp = await fetch(API_BASE + '/feedback', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          sequence: currentSequence,
          predicted_score: currentPrediction.suitability_score || 0,
          user_outcome: outcome,
          model_version: currentPrediction.model_version || null,
          email_for_followup: email || null,
          note: note || null,
        }),
      });
      if (!resp.ok) {
        // Re-enable buttons on failure
        feedbackBtns.forEach(b => b.disabled = false);
        return;
      }
      feedbackCard.hidden = true;
      feedbackThanks.hidden = false;
    } catch (err) {
      feedbackBtns.forEach(b => b.disabled = false);
    }
  });
});

sampleBtns.forEach(btn => {
  btn.addEventListener('click', () => {
    const key = btn.dataset.sample;
    const sample = SAMPLES[key];
    if (!sample) return;
    textarea.value = sample.seq;
    clearError();
    textarea.focus();
    textarea.scrollTop = 0;
  });
});

resetBtn.addEventListener('click', () => {
  resultsBox.hidden = true;
  clearError();
  textarea.value = '';
  textarea.focus();
  document.getElementById('predict').scrollIntoView({ behavior: 'smooth', block: 'start' });
});

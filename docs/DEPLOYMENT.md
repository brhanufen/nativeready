# NativeReady — Deployment Guide

Step-by-step instructions to take NativeReady from local code to a
publicly-accessible website. Total time: ~2 hours. Total cost: $5–25/month.

---

## Prerequisites

- The repo at `/Users/bfentaw2/startup/nativeready/` is in working order
- Local tests pass (`cd tests && python3 test_integration.py`)
- You have a domain name (or are willing to use a free `*.vercel.app` or `*.netlify.app` subdomain)

## Architecture (production)

```
            User browser
                │
                ▼
    nativeready.app  (Vercel / Netlify — free)
        ↑
        │ JS calls /predict
        ▼
    api.nativeready.app  (Railway / Render — $5/mo)
        │
        ▼
    model.joblib  (loaded into memory at boot)
```

---

## Step 1 — Push to GitHub (15 min)

GitHub becomes the source of truth and enables auto-deploys.

```bash
cd /Users/bfentaw2/startup/nativeready
git init
git add .
git commit -m "Initial NativeReady commit"
```

Create a private repo on github.com (e.g., `nativeready`), then:

```bash
git remote add origin git@github.com:YOURUSER/nativeready.git
git branch -M main
git push -u origin main
```

**Mitigation:** If you don't have GitHub CLI / SSH set up, use the GitHub
Desktop app instead. Same result, click-driven.

---

## Step 2 — Deploy backend to Railway (30 min)

Railway hosts the FastAPI service. The model file is bundled into the Docker image.

1. Sign up at https://railway.app (free trial; $5/month after)
2. Click **New Project** → **Deploy from GitHub repo** → select your repo
3. Set **Root Directory** to `backend`
4. Railway auto-detects the Dockerfile and builds
5. Add environment variables (none required for v1, but if you add Stripe later, this is where keys live)
6. Wait ~3 minutes for build + deploy
7. Click **Settings** → **Networking** → **Generate Domain** → you get something like `nativeready-api.up.railway.app`
8. Test it: `curl https://nativeready-api.up.railway.app/` should return `{"status":"ok","service":"nativeready"}`

**Mitigation if build fails:**
- Check that the `model/model.joblib` file is committed (run `git status` to verify)
- The Dockerfile assumes the model file is in `/app/model/model.joblib` — adjust path in `predictor.py` MODEL_PATH if needed for the deployed environment

**Mitigation for cold starts:**
- Free tier sleeps after inactivity. Pay $5/month for always-on.

---

## Step 3 — Deploy frontend to Netlify (15 min)

Netlify hosts the static frontend, free.

**Easiest path (drag-and-drop):**

1. Open `/Users/bfentaw2/startup/nativeready/frontend/script.js`
2. Update the `API_BASE` constant to point to your Railway URL:
   ```js
   const API_BASE = window.location.hostname === 'localhost'
     ? 'http://localhost:8000'
     : 'https://nativeready-api.up.railway.app';  // your Railway URL
   ```
3. Save
4. Go to https://app.netlify.com/drop
5. Drag the entire `frontend/` folder onto the page
6. Netlify gives you a URL like `https://wonderful-cat-12345.netlify.app`
7. Sign up to claim it permanently
8. Rename to `nativeready.netlify.app` (or similar) in Site settings

**Or via GitHub auto-deploy:**

1. Connect Netlify to your GitHub repo
2. Set **Base directory**: `frontend`
3. Set **Publish directory**: `frontend`
4. No build command needed
5. Auto-deploys on every git push

**Mitigation if API calls fail with CORS errors:**
- Verify that the backend has CORS enabled for your frontend origin
- The backend currently allows all origins. To restrict in production,
  edit `backend/main.py` and replace `allow_origins=["*"]` with your
  exact frontend URL.

---

## Step 4 — Custom domain (30 min, optional)

If you want `nativeready.app` instead of `nativeready.netlify.app`:

1. Buy domain at Cloudflare Registrar (cheapest, ~$15/year for `.app`)
2. In Netlify: Site settings → Domain management → Add custom domain
3. Netlify shows DNS records to add
4. In Cloudflare: DNS → Add the records Netlify gave you
5. Wait 5–30 minutes for DNS propagation
6. SSL cert is automatic via Netlify

**For the API subdomain (`api.nativeready.app`):**
1. In Railway: Settings → Networking → Custom Domain → enter `api.nativeready.app`
2. Railway gives you a CNAME target
3. In Cloudflare: add a CNAME record pointing `api` to that target
4. Update `frontend/script.js` `API_BASE` to `https://api.nativeready.app`
5. Redeploy frontend

**Mitigation: don't buy custom domain immediately.** Use the free
`*.netlify.app` URL until you've validated that anyone visits the site.
Saves $15/year if it doesn't work out.

---

## Step 5 — Monitoring (15 min)

Free monitoring so you know when things break.

**Sentry (free tier):**
1. Sign up at https://sentry.io
2. Create a Python project for the backend
3. Add the Sentry SDK to `backend/requirements.txt` and initialize in `main.py`
4. Errors now show up in the Sentry dashboard with stack traces

**Plausible or Vercel Analytics (free):**
- For the frontend, add a single `<script>` tag to `index.html`
- Tracks visitors and page views without cookies

**Uptime check (free):**
- https://uptimerobot.com — pings your `/` endpoint every 5 minutes
- Emails you if it's down

---

## Step 6 — Soft launch (1 weekend)

Once everything works at the public URL:

1. **Test it yourself once more** — visit the live site, paste a real sequence, verify the prediction looks right
2. **DM 5 specific people** in the native MS community: "I built this thing, would love your honest 2-minute review"
3. **Post on Twitter/X**, LinkedIn, r/Biochemistry — link the live site
4. **Monitor for 7 days** — answer every comment within 24 hours
5. **Track**: visitors, predictions made, errors, user feedback

---

## Step 7 — Iterate (ongoing)

Weekly cycle:

- Review feedback log → pick top 1–2 improvements
- Ship improvements via git push (auto-deploys on Vercel and Railway)
- Retrain model quarterly as new published native MS data accumulates
- Decision gates at day 30, 90, 180 (see `build_plan.docx` for details)

---

## Cost breakdown (first 6 months)

| Service       | Cost           | Required? |
|---------------|----------------|-----------|
| GitHub        | Free           | Yes       |
| Netlify       | Free           | Yes       |
| Railway       | $5/month       | Yes       |
| Cloudflare    | Free           | Yes       |
| Domain        | $15/year       | Optional  |
| Sentry        | Free tier      | Recommended |
| **Total**     | **~$45 + $15 domain** | |

---

## Troubleshooting

**Backend won't boot:**
- Check Railway logs: `railway logs` (or in the dashboard)
- Most common: missing dependency in `requirements.txt`

**Frontend can't reach backend:**
- Open browser dev tools → Network tab → look at the failed request
- Check `API_BASE` in `script.js` matches the deployed Railway URL
- CORS errors → backend needs to allow the frontend origin

**Model predictions look wrong:**
- Verify `model.joblib` was uploaded to Railway (check the build log)
- Verify `feature_names.json` is in `/app/model/` on Railway
- The backend's `extract_features()` MUST match the training pipeline exactly

**Costs higher than expected:**
- Set Railway billing alert at $20
- Free Netlify is unlimited bandwidth for static sites
- If model inference becomes the cost bottleneck, batch predictions or cache common sequences

---

*Last updated April 2026.*

# NativeReady — frontend

A static, single-page site for the NativeReady prediction tool. Three files do the
whole job: `index.html`, `style.css`, `script.js`. No build step. Host it anywhere
that can serve static files.

## Configure the API URL

Open `script.js` and edit the `API_BASE` constant near the top:

```js
const API_BASE = window.location.hostname === 'localhost'
  ? 'http://localhost:8000'
  : 'https://api.nativeready.app'; // change to your real backend
```

The frontend POSTs JSON `{ "sequence": "..." }` to `<API_BASE>/predict` and expects
the response shape documented in the project spec.

## Run locally

You can just open `index.html` in a browser, but the API call may need a real
origin. The simplest local server:

```bash
cd frontend
python3 -m http.server 5173
# then visit http://localhost:5173
```

When the page is served from `localhost`, `script.js` automatically points the
API at `http://localhost:8000` (run your FastAPI / Flask backend there).

## Deploy

Pick whichever is easiest.

### Netlify Drop (zero config)

1. Go to https://app.netlify.com/drop
2. Drag the `frontend/` folder onto the page.
3. Done. You get a URL like `https://nativeready-xyz.netlify.app`.

### GitHub Pages

1. Create a new GitHub repo, push the contents of `frontend/` to the root.
2. In the repo settings, enable Pages from the `main` branch, root folder.
3. Your site will be at `https://<your-username>.github.io/<repo-name>/`.

### Vercel

1. `npm i -g vercel` (one-time).
2. From inside `frontend/`, run `vercel` and accept the defaults.
3. Vercel detects a static site and deploys it.

### S3 + CloudFront

1. Create an S3 bucket configured for static website hosting.
2. Upload `index.html`, `style.css`, `script.js`.
3. Front it with CloudFront for HTTPS and a custom domain.

## CORS

Wherever you host the backend, make sure it sends
`Access-Control-Allow-Origin` for your frontend's origin. Otherwise the browser
will block the `/predict` call.

## What the page does

- Hero with the value prop.
- A textarea for FASTA / plain sequence input, with client-side cleaning
  (strips FASTA headers, whitespace, validates residue alphabet).
- Calls `POST /predict`, shows a loading spinner while waiting.
- Renders the suitability score, a colour-coded progress bar, risk-factor cards,
  and a recommendations list.
- Friendly error box if the API returns 4xx / 5xx or is unreachable.
- One-click example sequence (human carbonic anhydrase 2, UniProt P00918).
- Mobile-responsive down to ~375 px.

## Files

- `index.html` — markup and content
- `style.css`  — dark theme, cyan accent, matches the Traversa visual language
- `script.js`  — form handling, fetch, rendering, error UX

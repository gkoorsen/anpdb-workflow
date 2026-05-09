# Repository Bootstrap

## Local Git Setup

This project should be pushed as a clean source repository, not as a dump of all local data and generated outputs.

Recommended first commit contents:

- `.gitignore`
- `README.md`
- `docs/`
- `scripts/`
- `ppt/package.json`, `ppt/package-lock.json`, `ppt/build_ppt.js`, `ppt/build_ppt_old.js`

Do not commit:

- raw `data/` files
- `output/` contents
- MD trajectories/checkpoints
- `ppt/node_modules/`

## GitHub Creation

Create a new empty GitHub repository without adding a README, license, or `.gitignore` in the GitHub UI. Then run:

```bash
git remote add origin git@github.com:<user-or-org>/<repo>.git
git branch -M main
git add .gitignore README.md docs scripts ppt/package.json ppt/package-lock.json ppt/build_ppt.js ppt/build_ppt_old.js
git commit -m "Prepare repository for GPU MD workflow"
git push -u origin main
```

If HTTPS is preferred:

```bash
git remote add origin https://github.com/<user-or-org>/<repo>.git
```

## WSL GPU Machine Setup

On the NVIDIA GPU machine:

```bash
git clone git@github.com:<user-or-org>/<repo>.git
cd <repo>
```

The production MD workflow should be added as source-controlled scripts/configs, while large inputs and outputs should be staged outside Git, for example:

```bash
mkdir -p data md_runs
```

The next implementation step is to add a dedicated GPU MD environment file and production runner with CUDA checks, checkpointing, resumable 100 ns runs, and run manifests.

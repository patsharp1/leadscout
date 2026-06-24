#!/usr/bin/env bash
# setup_github.sh — initialise the Leadscout repo and push to GitHub Pages
# Usage: bash setup_github.sh <github-username>
# Example: bash setup_github.sh patsharp

set -e

OUTPUTS_DIR="$(cd "$(dirname "$0")" && pwd)"
USERNAME="${1:-}"

if [ -z "$USERNAME" ]; then
  echo "Usage: bash setup_github.sh <your-github-username>"
  exit 1
fi

REPO_NAME="leadscout"
REMOTE_URL="https://github.com/${USERNAME}/${REPO_NAME}.git"

echo ""
echo "═══════════════════════════════════════════"
echo "  Leadscout → GitHub Pages setup"
echo "  Repo: ${REMOTE_URL}"
echo "═══════════════════════════════════════════"
echo ""

cd "$OUTPUTS_DIR"

# Initialise git if not already done
if [ ! -d ".git" ]; then
  git init
  git checkout -b main
  echo "✓ git init"
else
  echo "✓ existing git repo found"
  git checkout -B main 2>/dev/null || true
fi

# Stage everything (respects .gitignore)
git add .
echo "✓ files staged"

# Commit
git commit -m "feat: initial Leadscout — competitive chemistry intelligence feed" \
  --allow-empty
echo "✓ committed"

# Set remote (replace if exists)
git remote remove origin 2>/dev/null || true
git remote add origin "$REMOTE_URL"
echo "✓ remote set to $REMOTE_URL"

# Push
echo ""
echo "Pushing to GitHub (you may be prompted for credentials)…"
git push -u origin main --force
echo ""
echo "═══════════════════════════════════════════"
echo "  ✓ Code pushed!"
echo ""
echo "  NEXT: enable GitHub Pages"
echo "  1. Go to: https://github.com/${USERNAME}/${REPO_NAME}/settings/pages"
echo "  2. Source → 'Deploy from a branch'"
echo "  3. Branch → main  /  (root)"
echo "  4. Click Save"
echo ""
echo "  Your site will be live at:"
echo "  https://${USERNAME}.github.io/${REPO_NAME}"
echo "═══════════════════════════════════════════"

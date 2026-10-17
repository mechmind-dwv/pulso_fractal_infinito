#!/data/data/com.termux/files/usr/bin/bash
set -e
cd "$(dirname "$0")/.."

if [ ! -d .git ]; then
    echo "No hay repositorio git aca. Inicializando..."
    git init
    git branch -M main
fi

if ! git remote | grep -q origin; then
    echo "No hay remoto 'origin' configurado. Corre esto una vez:"
    echo "  git remote add origin https://github.com/TU-ORG/TU-REPO.git"
    exit 1
fi

MSG="${1:-actualizacion $(date '+%Y-%m-%d %H:%M')}"
git add .
git commit -m "$MSG" || echo "(nada nuevo que commitear)"
git push origin main

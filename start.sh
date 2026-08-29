#!/usr/bin/env bash
# Serve the tool and open it in your browser.
set -euo pipefail
cd "$(dirname "$0")"
PORT="${PORT:-8123}"
PY=$(command -v python3 || command -v python)

echo "Warehouse Area Measure -> http://localhost:$PORT"
echo "Press Ctrl-C to stop."

"$PY" -m http.server "$PORT" &
SERVER=$!
trap 'kill $SERVER 2>/dev/null || true' EXIT
sleep 1

if   command -v xdg-open >/dev/null; then xdg-open "http://localhost:$PORT"
elif command -v open     >/dev/null; then open     "http://localhost:$PORT"
fi
wait $SERVER

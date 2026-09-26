@echo off
rem Launch Edge with CDP debug port using the DEFAULT profile,
rem so your existing GitHub login is preserved.
rem IMPORTANT: close ALL Edge windows first, otherwise Edge ignores this flag.
start "" "C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe" --remote-debugging-port=9222 --remote-allow-origins=* --restore-last-session --no-first-run --no-default-browser-check "https://github.com"
echo.
echo Edge launched with debug port 9222 (default profile, login preserved).
echo If port 9222 does not respond, close every Edge window and run this again.

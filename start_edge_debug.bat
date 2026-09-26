@echo off
rem Launch Edge with a CDP debug port so the assistant can drive it.
rem Uses a separate profile folder; your normal Edge profile is not touched.
start "" "C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe" --remote-debugging-port=9222 --user-data-dir="D:\edge-debug-profile" --no-first-run --no-default-browser-check "https://github.com"
echo Edge launched with debug port 9222.

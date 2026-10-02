#!/usr/bin/env python3
"""Embed data/xy_data.json into viewer_template.html -> index.html."""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(HERE, "viewer_template.html")) as f:
    html = f.read()
with open(os.path.join(HERE, "data", "xy_data.json")) as f:
    data = f.read()
with open(os.path.join(HERE, "index.html"), "w") as f:
    f.write(html.replace("/*DATA*/null", data.replace("</", "<\\/")))
print("wrote index.html")

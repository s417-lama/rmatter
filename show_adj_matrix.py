import os
import sys
import webbrowser
import numpy as np
import plotly.graph_objects as go
import plotly.express as px

if len(sys.argv) < 2:
    print("Usage: python3 show_adj_matrix.py <input_file> (<max_cells>)")
    sys.exit(1)

input_file = sys.argv[1]

max_cells = int(sys.argv[2]) if len(sys.argv) >= 3 else 1024

edges = np.loadtxt(input_file, dtype=int)

n = np.max(edges) + 1

n_cells = min(n, max_cells)

heatmap_data, _, _ = np.histogram2d(edges[:, 0], edges[:, 1], bins=n_cells, range=[[0, n-1], [0, n-1]])

colorscheme = px.colors.sequential.Magma
colorscale = [[4**(-(len(colorscheme) - i - 1)), c] for i, c in enumerate(colorscheme)]
colorscale[0] = [0, colorscheme[0]]

fig = go.Figure(data=go.Heatmap(
    z=heatmap_data,
    colorscale=colorscale,
    showscale=False,
))

fig.update_layout(
    width=500,
    height=500,
    margin=dict(t=50, b=50, l=50, r=50),
    xaxis=dict(showline=True, linecolor="grey", mirror=True, range=[0, n_cells-1], tickvals=[0, n_cells-1], ticktext=[0, n], side="top"),
    yaxis=dict(showline=True, linecolor="grey", mirror=True, range=[n_cells-1, 0], tickvals=[0, n_cells-1], ticktext=[0, n]),
)

output_file = input_file + ".html"
fig.write_html(output_file, include_plotlyjs="cdn", config={"toImageButtonOptions" : {"format" : "svg"}})
webbrowser.open("file://" + os.path.realpath(output_file))

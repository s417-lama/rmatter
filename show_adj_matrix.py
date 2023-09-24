import os
import sys
import webbrowser
import dask.array as da
import dask.dataframe as dd
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio

if len(sys.argv) < 2:
    print("Usage: python3 show_adj_matrix.py <input_file> (<max_hist2d_bins>) (<max_hist1d_bins>)")
    sys.exit(1)

input_file = sys.argv[1]
max_bins_2d = int(sys.argv[2]) if len(sys.argv) >= 3 else 1024
max_bins_1d = int(sys.argv[3]) if len(sys.argv) >= 4 else 128

edges = dd.read_csv(input_file, header=None, sep="\s+", dtype=int).values

n = edges.max().compute() + 1
n_bins_2d = min(n, max_bins_2d)
n_bins_1d = min(n_bins_2d, max_bins_1d)

hist2d, _, _ = da.histogram2d(edges[:, 0], edges[:, 1], bins=n_bins_2d, range=((0, n-1), (0, n-1)))
hist1d_src, _ = da.histogram(edges[:, 0], bins=n_bins_1d, range=(0, n-1))
hist1d_dst, _ = da.histogram(edges[:, 1], bins=n_bins_1d, range=(0, n-1))

colorscheme = px.colors.sequential.Magma
colorscale = [[4**(-(len(colorscheme) - i - 1)), c] for i, c in enumerate(colorscheme)]
colorscale[0] = [0, colorscheme[0]]

pio.templates.default = "plotly_white"

fig = go.Figure().set_subplots(
    rows=2, cols=2,
    column_widths=[0.82, 0.18],
    row_heights=[0.82, 0.18],
    horizontal_spacing = 0.08,
    vertical_spacing = 0.08,
)

fig.add_trace(go.Heatmap(
    z=hist2d,
    colorscale=colorscale,
    showscale=False,
), row=1, col=1)

fig.add_trace(go.Bar(
    x=hist1d_src,
    orientation="h",
    marker=dict(color=colorscheme[len(colorscheme) // 2]),
    width=1,
), row=1, col=2)

fig.add_trace(go.Bar(
    y=hist1d_dst,
    marker=dict(color=colorscheme[len(colorscheme) // 2]),
    width=1,
), row=2, col=1)

fig.update_layout(
    width=500,
    height=500,
    margin=dict(t=50, b=50, l=50, r=50),
    xaxis=dict(showline=True, linecolor="grey", mirror=True, range=[-0.5, n_bins_2d-0.5], tickvals=[0, n_bins_2d-1], ticktext=[0, n], side="top"),
    yaxis=dict(showline=True, linecolor="grey", mirror=True, range=[n_bins_2d-0.5, -0.5], tickvals=[0, n_bins_2d-1], ticktext=[0, n]),
    xaxis2=dict(showline=True, linecolor="grey", mirror=True, showticklabels=True, ticks="outside", side="top"),
    yaxis2=dict(showline=True, linecolor="grey", mirror=True, showticklabels=False, range=[n_bins_1d-0.5, -0.5]),
    xaxis3=dict(showline=True, linecolor="grey", mirror=True, showticklabels=False, range=[-0.5, n_bins_1d-0.5]),
    yaxis3=dict(showline=True, linecolor="grey", mirror=True, showticklabels=True, ticks="outside", autorange="reversed"),
    showlegend=False,
)

output_file = input_file + ".html"
fig.write_html(output_file, include_plotlyjs="cdn", config={"toImageButtonOptions" : {"format" : "svg"}})
webbrowser.open("file://" + os.path.realpath(output_file))

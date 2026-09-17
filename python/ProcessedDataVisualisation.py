import numpy as np
import pandas as pd
from dash import Dash, dcc, html, Input, Output, ALL
import plotly.express as px
import dash_bootstrap_components as dbc


class ProcessedDataVisualisation():
    def __init__(self, fname):
        self.data = pd.read_csv(fname)
        self.app = Dash(
            external_stylesheets=[dbc.themes.SIMPLEX],
            suppress_callback_exceptions=True
        )
        self.build_layout()
        self.add_callbacks()

    def build_layout(self):
        self.numeric_cols = self.data.select_dtypes(
            include="number").columns.tolist()
        self.numeric_cols = [
            s for s in self.numeric_cols if not s.startswith("Unnamed")]
        print(self.numeric_cols)
        self.categorical_cols = self.data.select_dtypes(
            exclude="number").columns.tolist()
        self.all_cols = self.data.columns.tolist()
        self.all_cols = [
            s for s in self.all_cols if not s.startswith("Unnamed")]

        self.app.title = "Scatter Explorer"

        def axis_dropdown(id_, default, options=self.all_cols, clearable=False):
            return dcc.Dropdown(
                id=id_,
                options=[{"label": c, "value": c} for c in options],
                value=default,
                clearable=clearable,
                style={"marginBottom": "10px"},
            )

        def make_filter_control(col):
            """Build a RangeSlider for numeric columns or a Checklist for categoricals.
            Every control uses a pattern-matching id so one callback can read all of them."""
            if col in self.numeric_cols:
                lo, hi = float(self.data[col].min()), float(
                    self.data[col].max())
                step = (hi - lo) / 100 if hi > lo else 1
                return html.Div(
                    [
                        html.Label(col, style={"fontWeight": "bold"}),
                        dcc.RangeSlider(
                            id={"type": "filter-numeric", "col": col},
                            min=lo,
                            max=hi,
                            step=step,
                            value=[lo, hi],
                            tooltip={"placement": "bottom",
                                     "always_visible": False},
                            allowCross=False,
                        ),
                    ],
                    style={"marginBottom": "24px"},
                )
            else:
                opts = sorted(
                    self.data[col].dropna().unique().tolist(), key=str)
                return html.Div(
                    [
                        html.Label(col, style={"fontWeight": "bold"}),
                        dcc.Checklist(
                            id={"type": "filter-categorical", "col": col},
                            options=[{"label": str(o), "value": o}
                                     for o in opts],
                            value=opts,
                            inline=True,
                            style={"marginTop": "4px"},
                        ),
                    ],
                    style={"marginBottom": "24px"},
                )

        self.filter_controls = [make_filter_control(c) for c in self.all_cols]

        self.app.layout = html.Div(
            style={"display": "flex", "fontFamily": "Arial, sans-serif"},
            children=[
                # --- Sidebar controls ---------------------------------------------
                html.Div(
                    style={
                        "width": "320px",
                        "padding": "20px",
                        "borderRight": "1px solid #ddd",
                        "height": "100vh",
                        "overflowY": "auto",
                        "boxSizing": "border-box",
                    },
                    children=[
                        html.H3("Axes & Style"),
                        html.Label("X axis"),
                        axis_dropdown(
                            "x-axis", self.numeric_cols[0] if self.numeric_cols else self.all_cols[0]),
                        html.Label("Y axis"),
                        axis_dropdown(
                            "y-axis",
                            self.numeric_cols[1] if len(
                                self.numeric_cols) > 1 else self.numeric_cols[0],
                        ),
                        html.Label("Color"),
                        axis_dropdown(
                            "color-axis",
                            self.categorical_cols[0] if self.categorical_cols else None,
                            options=self.all_cols,
                            clearable=True,
                        ),
                        html.Label("Size"),
                        axis_dropdown(
                            "size-axis",
                            self.numeric_cols[-1] if self.numeric_cols else None,
                            options=self.numeric_cols,
                            clearable=True,
                        ),
                        html.Hr(),
                        html.H3("Filters"),
                        html.Div(id="filters-container",
                                 children=self.filter_controls),
                    ],
                ),
                # --- Plot ------------------------------------------------------------
                html.Div(
                    style={"flex": "1", "padding": "20px",
                           "boxSizing": "border-box"},
                    children=[
                        html.Div(id="row-count",
                                 style={"marginBottom": "10px", "color": "#555"}),
                        dcc.Graph(id="scatter-plot", style={"height": "85vh"}),
                    ],
                ),
            ],
        )

    def add_callbacks(self):
        @self.app.callback(
            Output("scatter-plot", "figure"),
            Output("row-count", "children"),
            Input("x-axis", "value"),
            Input("y-axis", "value"),
            Input("color-axis", "value"),
            Input("size-axis", "value"),
            Input({"type": "filter-numeric", "col": ALL}, "value"),
            Input({"type": "filter-categorical", "col": ALL}, "value"),
        )
        def update_plot(x_col, y_col, color_col, size_col, numeric_values, categorical_values):
            dff = self.data.copy()

            # Apply numeric range filters (order matches the order numeric columns appear in all_cols).
            numeric_filter_cols = [
                c for c in self.all_cols if c in self.numeric_cols]
            for col, (lo, hi) in zip(numeric_filter_cols, numeric_values):
                dff = dff[(dff[col] >= lo) & (dff[col] <= hi)]

            # Apply categorical filters.
            categorical_filter_cols = [
                c for c in self.all_cols if c in self.categorical_cols]
            for col, selected in zip(categorical_filter_cols, categorical_values):
                dff = dff[dff[col].isin(selected)]

            scatter_kwargs = dict(x=x_col, y=y_col)
            if color_col:
                scatter_kwargs["color"] = color_col
            if size_col:
                # Plotly's size scale requires non-negative values.
                dff = dff[dff[size_col].notna()]
                dff = dff[dff[size_col] >= 0]
                scatter_kwargs["size"] = size_col

            if dff.empty:
                fig = px.scatter(title="No data matches the current filters")
            else:
                fig = px.scatter(dff, **scatter_kwargs, opacity=0.75)
                fig.update_layout(template="plotly_white",
                                  margin=dict(l=40, r=20, t=40, b=40))

            count_text = f"Showing {len(dff)} of {len(self.data)} rows"
            return fig, count_text

    def run(self):
        self.app.run(debug=True)


if __name__ == "__main__":
    visuals = ProcessedDataVisualisation(
        "./results_voice_processed_data/processed_data.csv")
    visuals.run()

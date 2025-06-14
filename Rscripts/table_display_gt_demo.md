table_display_gt_demo
================
Janet Young

2025-05-27

``` r
# Define the start and end dates for the data range
start_date <- "2010-06-07"
end_date <- "2010-06-14"

# Create a gt table based on preprocessed
# `sp500` table data
sp500 |>
    dplyr::filter(date >= start_date & date <= end_date) |>
    dplyr::select(-adj_close) |>
    gt() |>
    tab_header(
        title = "S&P 500",
        subtitle = glue::glue("{start_date} to {end_date}")
    ) |>
    fmt_currency() |>
    fmt_date(columns = date, date_style = "wd_m_day_year") |>
    fmt_number(columns = volume, suffixing = TRUE)
```

<div id="cisqbzieai" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#cisqbzieai table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#cisqbzieai thead, #cisqbzieai tbody, #cisqbzieai tfoot, #cisqbzieai tr, #cisqbzieai td, #cisqbzieai th {
  border-style: none;
}
&#10;#cisqbzieai p {
  margin: 0;
  padding: 0;
}
&#10;#cisqbzieai .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}
&#10;#cisqbzieai .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#cisqbzieai .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}
&#10;#cisqbzieai .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}
&#10;#cisqbzieai .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#cisqbzieai .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#cisqbzieai .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#cisqbzieai .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}
&#10;#cisqbzieai .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}
&#10;#cisqbzieai .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#cisqbzieai .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#cisqbzieai .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}
&#10;#cisqbzieai .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#cisqbzieai .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}
&#10;#cisqbzieai .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}
&#10;#cisqbzieai .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#cisqbzieai .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#cisqbzieai .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}
&#10;#cisqbzieai .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#cisqbzieai .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}
&#10;#cisqbzieai .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#cisqbzieai .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#cisqbzieai .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#cisqbzieai .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#cisqbzieai .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#cisqbzieai .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#cisqbzieai .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#cisqbzieai .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#cisqbzieai .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#cisqbzieai .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#cisqbzieai .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#cisqbzieai .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#cisqbzieai .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#cisqbzieai .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#cisqbzieai .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#cisqbzieai .gt_left {
  text-align: left;
}
&#10;#cisqbzieai .gt_center {
  text-align: center;
}
&#10;#cisqbzieai .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#cisqbzieai .gt_font_normal {
  font-weight: normal;
}
&#10;#cisqbzieai .gt_font_bold {
  font-weight: bold;
}
&#10;#cisqbzieai .gt_font_italic {
  font-style: italic;
}
&#10;#cisqbzieai .gt_super {
  font-size: 65%;
}
&#10;#cisqbzieai .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#cisqbzieai .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#cisqbzieai .gt_indent_1 {
  text-indent: 5px;
}
&#10;#cisqbzieai .gt_indent_2 {
  text-indent: 10px;
}
&#10;#cisqbzieai .gt_indent_3 {
  text-indent: 15px;
}
&#10;#cisqbzieai .gt_indent_4 {
  text-indent: 20px;
}
&#10;#cisqbzieai .gt_indent_5 {
  text-indent: 25px;
}
&#10;#cisqbzieai .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}
&#10;#cisqbzieai div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_heading">
      <td colspan="6" class="gt_heading gt_title gt_font_normal" style>S&amp;P 500</td>
    </tr>
    <tr class="gt_heading">
      <td colspan="6" class="gt_heading gt_subtitle gt_font_normal gt_bottom_border" style>2010-06-07 to 2010-06-14</td>
    </tr>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="date">date</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="open">open</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="high">high</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="low">low</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="close">close</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="volume">volume</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="date" class="gt_row gt_right">Mon, Jun 14, 2010</td>
<td headers="open" class="gt_row gt_right">$1,095.00</td>
<td headers="high" class="gt_row gt_right">$1,105.91</td>
<td headers="low" class="gt_row gt_right">$1,089.03</td>
<td headers="close" class="gt_row gt_right">$1,089.63</td>
<td headers="volume" class="gt_row gt_right">4.43B</td></tr>
    <tr><td headers="date" class="gt_row gt_right">Fri, Jun 11, 2010</td>
<td headers="open" class="gt_row gt_right">$1,082.65</td>
<td headers="high" class="gt_row gt_right">$1,092.25</td>
<td headers="low" class="gt_row gt_right">$1,077.12</td>
<td headers="close" class="gt_row gt_right">$1,091.60</td>
<td headers="volume" class="gt_row gt_right">4.06B</td></tr>
    <tr><td headers="date" class="gt_row gt_right">Thu, Jun 10, 2010</td>
<td headers="open" class="gt_row gt_right">$1,058.77</td>
<td headers="high" class="gt_row gt_right">$1,087.85</td>
<td headers="low" class="gt_row gt_right">$1,058.77</td>
<td headers="close" class="gt_row gt_right">$1,086.84</td>
<td headers="volume" class="gt_row gt_right">5.14B</td></tr>
    <tr><td headers="date" class="gt_row gt_right">Wed, Jun 9, 2010</td>
<td headers="open" class="gt_row gt_right">$1,062.75</td>
<td headers="high" class="gt_row gt_right">$1,077.74</td>
<td headers="low" class="gt_row gt_right">$1,052.25</td>
<td headers="close" class="gt_row gt_right">$1,055.69</td>
<td headers="volume" class="gt_row gt_right">5.98B</td></tr>
    <tr><td headers="date" class="gt_row gt_right">Tue, Jun 8, 2010</td>
<td headers="open" class="gt_row gt_right">$1,050.81</td>
<td headers="high" class="gt_row gt_right">$1,063.15</td>
<td headers="low" class="gt_row gt_right">$1,042.17</td>
<td headers="close" class="gt_row gt_right">$1,062.00</td>
<td headers="volume" class="gt_row gt_right">6.19B</td></tr>
    <tr><td headers="date" class="gt_row gt_right">Mon, Jun 7, 2010</td>
<td headers="open" class="gt_row gt_right">$1,065.84</td>
<td headers="high" class="gt_row gt_right">$1,071.36</td>
<td headers="low" class="gt_row gt_right">$1,049.86</td>
<td headers="close" class="gt_row gt_right">$1,050.47</td>
<td headers="volume" class="gt_row gt_right">5.47B</td></tr>
  </tbody>
  &#10;  
</table>
</div>

Show sequence alignment as a table - play around with this as an
(alternative to ggmsa:

``` r
## this is actually a slice of Orf9b prot seq from SARS-CoV2 and other viruses
seq_slice <- c("FQLT",
               "FQLI",
               "FQLT",
               "FRLT",
               "FRLI",
               "FRLT",
               "FRLT",
               "FQLT",
               "FQST",
               "FQST")
seq_slice_tbl <- seq_slice %>% 
    strsplit("") %>% 
    as.data.frame() %>% 
    set_names(nm=paste0("seq", 1:10)) %>%
    t() %>% 
    as.data.frame() %>%
    set_names(nm=paste0("pos", 1:4)) %>%
    as_tibble(rownames="id")
```

``` r
seq_slice_tbl %>% 
    ## gt alone makes a decent-looking table
    gt() %>% 
    ## a bunch of formatting things:
    cols_align(align = "center", columns=-id) %>% 
    opt_table_font(font = list(google_font(name = "Courier"))) %>% 
    rm_header() %>% 
    tab_options(column_labels.hidden = TRUE,
                table_body.border.bottom.style = "hidden",
                table_body.hlines.style = "hidden",
                data_row.padding=0,
                data_row.padding.horizontal=2,
                ## row striping not visible until I knit. Options don't seem to work. Might be an Rstudio bug
                row.striping.include_table_body = FALSE,
                row.striping.include_stub=FALSE) %>% 
    ## px is number of pixels
    cols_width(id ~ px(75)) 
```

<div id="aljqzqidae" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>@import url("https://fonts.googleapis.com/css2?family=Courier:ital,wght@0,100;0,200;0,300;0,400;0,500;0,600;0,700;0,800;0,900;1,100;1,200;1,300;1,400;1,500;1,600;1,700;1,800;1,900&display=swap");
#aljqzqidae table {
  font-family: Courier, system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#aljqzqidae thead, #aljqzqidae tbody, #aljqzqidae tfoot, #aljqzqidae tr, #aljqzqidae td, #aljqzqidae th {
  border-style: none;
}
&#10;#aljqzqidae p {
  margin: 0;
  padding: 0;
}
&#10;#aljqzqidae .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}
&#10;#aljqzqidae .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#aljqzqidae .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}
&#10;#aljqzqidae .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}
&#10;#aljqzqidae .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#aljqzqidae .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#aljqzqidae .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#aljqzqidae .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}
&#10;#aljqzqidae .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}
&#10;#aljqzqidae .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#aljqzqidae .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#aljqzqidae .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}
&#10;#aljqzqidae .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#aljqzqidae .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}
&#10;#aljqzqidae .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}
&#10;#aljqzqidae .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#aljqzqidae .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#aljqzqidae .gt_row {
  padding-top: 0px;
  padding-bottom: 0px;
  padding-left: 2px;
  padding-right: 2px;
  margin: 10px;
  border-top-style: hidden;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}
&#10;#aljqzqidae .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 2px;
  padding-right: 2px;
}
&#10;#aljqzqidae .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}
&#10;#aljqzqidae .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#aljqzqidae .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#aljqzqidae .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#aljqzqidae .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#aljqzqidae .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#aljqzqidae .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#aljqzqidae .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#aljqzqidae .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#aljqzqidae .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#aljqzqidae .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#aljqzqidae .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: hidden;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#aljqzqidae .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#aljqzqidae .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#aljqzqidae .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#aljqzqidae .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#aljqzqidae .gt_left {
  text-align: left;
}
&#10;#aljqzqidae .gt_center {
  text-align: center;
}
&#10;#aljqzqidae .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#aljqzqidae .gt_font_normal {
  font-weight: normal;
}
&#10;#aljqzqidae .gt_font_bold {
  font-weight: bold;
}
&#10;#aljqzqidae .gt_font_italic {
  font-style: italic;
}
&#10;#aljqzqidae .gt_super {
  font-size: 65%;
}
&#10;#aljqzqidae .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#aljqzqidae .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#aljqzqidae .gt_indent_1 {
  text-indent: 5px;
}
&#10;#aljqzqidae .gt_indent_2 {
  text-indent: 10px;
}
&#10;#aljqzqidae .gt_indent_3 {
  text-indent: 15px;
}
&#10;#aljqzqidae .gt_indent_4 {
  text-indent: 20px;
}
&#10;#aljqzqidae .gt_indent_5 {
  text-indent: 25px;
}
&#10;#aljqzqidae .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}
&#10;#aljqzqidae div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" style="table-layout:fixed;" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <colgroup>
    <col style="width:75px;"/>
    <col/>
    <col/>
    <col/>
    <col/>
  </colgroup>
  &#10;  <tbody class="gt_table_body">
    <tr><td headers="id" class="gt_row gt_left">seq1</td>
<td headers="pos1" class="gt_row gt_center">F</td>
<td headers="pos2" class="gt_row gt_center">Q</td>
<td headers="pos3" class="gt_row gt_center">L</td>
<td headers="pos4" class="gt_row gt_center">T</td></tr>
    <tr><td headers="id" class="gt_row gt_left">seq2</td>
<td headers="pos1" class="gt_row gt_center">F</td>
<td headers="pos2" class="gt_row gt_center">Q</td>
<td headers="pos3" class="gt_row gt_center">L</td>
<td headers="pos4" class="gt_row gt_center">I</td></tr>
    <tr><td headers="id" class="gt_row gt_left">seq3</td>
<td headers="pos1" class="gt_row gt_center">F</td>
<td headers="pos2" class="gt_row gt_center">Q</td>
<td headers="pos3" class="gt_row gt_center">L</td>
<td headers="pos4" class="gt_row gt_center">T</td></tr>
    <tr><td headers="id" class="gt_row gt_left">seq4</td>
<td headers="pos1" class="gt_row gt_center">F</td>
<td headers="pos2" class="gt_row gt_center">R</td>
<td headers="pos3" class="gt_row gt_center">L</td>
<td headers="pos4" class="gt_row gt_center">T</td></tr>
    <tr><td headers="id" class="gt_row gt_left">seq5</td>
<td headers="pos1" class="gt_row gt_center">F</td>
<td headers="pos2" class="gt_row gt_center">R</td>
<td headers="pos3" class="gt_row gt_center">L</td>
<td headers="pos4" class="gt_row gt_center">I</td></tr>
    <tr><td headers="id" class="gt_row gt_left">seq6</td>
<td headers="pos1" class="gt_row gt_center">F</td>
<td headers="pos2" class="gt_row gt_center">R</td>
<td headers="pos3" class="gt_row gt_center">L</td>
<td headers="pos4" class="gt_row gt_center">T</td></tr>
    <tr><td headers="id" class="gt_row gt_left">seq7</td>
<td headers="pos1" class="gt_row gt_center">F</td>
<td headers="pos2" class="gt_row gt_center">R</td>
<td headers="pos3" class="gt_row gt_center">L</td>
<td headers="pos4" class="gt_row gt_center">T</td></tr>
    <tr><td headers="id" class="gt_row gt_left">seq8</td>
<td headers="pos1" class="gt_row gt_center">F</td>
<td headers="pos2" class="gt_row gt_center">Q</td>
<td headers="pos3" class="gt_row gt_center">L</td>
<td headers="pos4" class="gt_row gt_center">T</td></tr>
    <tr><td headers="id" class="gt_row gt_left">seq9</td>
<td headers="pos1" class="gt_row gt_center">F</td>
<td headers="pos2" class="gt_row gt_center">Q</td>
<td headers="pos3" class="gt_row gt_center">S</td>
<td headers="pos4" class="gt_row gt_center">T</td></tr>
    <tr><td headers="id" class="gt_row gt_left">seq10</td>
<td headers="pos1" class="gt_row gt_center">F</td>
<td headers="pos2" class="gt_row gt_center">Q</td>
<td headers="pos3" class="gt_row gt_center">S</td>
<td headers="pos4" class="gt_row gt_center">T</td></tr>
  </tbody>
  &#10;  
</table>
</div>

# Finished

``` r
sessionInfo()
```

    ## R version 4.4.2 (2024-10-31)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Sequoia 15.5
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/Los_Angeles
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] gt_1.0.0        here_1.0.1      lubridate_1.9.4 forcats_1.0.0  
    ##  [5] stringr_1.5.1   dplyr_1.1.4     purrr_1.0.2     readr_2.1.5    
    ##  [9] tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1   tidyverse_2.0.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.6      compiler_4.4.2    tidyselect_1.2.1  xml2_1.3.6       
    ##  [5] scales_1.3.0      yaml_2.3.10       fastmap_1.2.0     R6_2.5.1         
    ##  [9] generics_0.1.3    bigD_0.3.1        knitr_1.49        rprojroot_2.0.4  
    ## [13] munsell_0.5.1     pillar_1.10.1     tzdb_0.4.0        rlang_1.1.4      
    ## [17] stringi_1.8.4     xfun_0.50         sass_0.4.9        timechange_0.3.0 
    ## [21] cli_3.6.3         withr_3.0.2       magrittr_2.0.3    digest_0.6.37    
    ## [25] grid_4.4.2        rstudioapi_0.17.1 hms_1.1.3         lifecycle_1.0.4  
    ## [29] vctrs_0.6.5       evaluate_1.0.3    glue_1.8.0        colorspace_2.1-1 
    ## [33] rmarkdown_2.29    tools_4.4.2       pkgconfig_2.0.3   htmltools_0.5.8.1

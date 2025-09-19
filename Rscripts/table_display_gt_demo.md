table_display_gt_demo
================
Janet Young

2025-09-18

Show without adding as_raw_html()

``` r
# Define the start and end dates for the data range
start_date <- "2010-06-07"
end_date <- "2010-06-14"

# Create a gt table based on preprocessed `sp500` table data
sp500_gt <- sp500 |>
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
sp500_gt
```

<div id="jnaoqudxok" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#jnaoqudxok table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#jnaoqudxok thead, #jnaoqudxok tbody, #jnaoqudxok tfoot, #jnaoqudxok tr, #jnaoqudxok td, #jnaoqudxok th {
  border-style: none;
}
&#10;#jnaoqudxok p {
  margin: 0;
  padding: 0;
}
&#10;#jnaoqudxok .gt_table {
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
&#10;#jnaoqudxok .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#jnaoqudxok .gt_title {
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
&#10;#jnaoqudxok .gt_subtitle {
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
&#10;#jnaoqudxok .gt_heading {
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
&#10;#jnaoqudxok .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#jnaoqudxok .gt_col_headings {
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
&#10;#jnaoqudxok .gt_col_heading {
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
&#10;#jnaoqudxok .gt_column_spanner_outer {
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
&#10;#jnaoqudxok .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#jnaoqudxok .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#jnaoqudxok .gt_column_spanner {
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
&#10;#jnaoqudxok .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#jnaoqudxok .gt_group_heading {
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
&#10;#jnaoqudxok .gt_empty_group_heading {
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
&#10;#jnaoqudxok .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#jnaoqudxok .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#jnaoqudxok .gt_row {
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
&#10;#jnaoqudxok .gt_stub {
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
&#10;#jnaoqudxok .gt_stub_row_group {
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
&#10;#jnaoqudxok .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#jnaoqudxok .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#jnaoqudxok .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#jnaoqudxok .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#jnaoqudxok .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#jnaoqudxok .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#jnaoqudxok .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#jnaoqudxok .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#jnaoqudxok .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#jnaoqudxok .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#jnaoqudxok .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#jnaoqudxok .gt_footnotes {
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
&#10;#jnaoqudxok .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#jnaoqudxok .gt_sourcenotes {
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
&#10;#jnaoqudxok .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#jnaoqudxok .gt_left {
  text-align: left;
}
&#10;#jnaoqudxok .gt_center {
  text-align: center;
}
&#10;#jnaoqudxok .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#jnaoqudxok .gt_font_normal {
  font-weight: normal;
}
&#10;#jnaoqudxok .gt_font_bold {
  font-weight: bold;
}
&#10;#jnaoqudxok .gt_font_italic {
  font-style: italic;
}
&#10;#jnaoqudxok .gt_super {
  font-size: 65%;
}
&#10;#jnaoqudxok .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#jnaoqudxok .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#jnaoqudxok .gt_indent_1 {
  text-indent: 5px;
}
&#10;#jnaoqudxok .gt_indent_2 {
  text-indent: 10px;
}
&#10;#jnaoqudxok .gt_indent_3 {
  text-indent: 15px;
}
&#10;#jnaoqudxok .gt_indent_4 {
  text-indent: 20px;
}
&#10;#jnaoqudxok .gt_indent_5 {
  text-indent: 25px;
}
&#10;#jnaoqudxok .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}
&#10;#jnaoqudxok div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
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

Show including as_raw_html()

``` r
sp500_gt |>
    ## if I don't include as_raw_html(), when I render and display as a github page I get a lot of extra stuff above the table
    as_raw_html()
```

<div id="zrysstqbds" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
  &#10;  <table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false" style="-webkit-font-smoothing: antialiased; -moz-osx-font-smoothing: grayscale; font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji'; display: table; border-collapse: collapse; line-height: normal; margin-left: auto; margin-right: auto; color: #333333; font-size: 16px; font-weight: normal; font-style: normal; background-color: #FFFFFF; width: auto; border-top-style: solid; border-top-width: 2px; border-top-color: #A8A8A8; border-right-style: none; border-right-width: 2px; border-right-color: #D3D3D3; border-bottom-style: solid; border-bottom-width: 2px; border-bottom-color: #A8A8A8; border-left-style: none; border-left-width: 2px; border-left-color: #D3D3D3;" bgcolor="#FFFFFF">
  <thead style="border-style: none;">
    <tr class="gt_heading" style="border-style: none; background-color: #FFFFFF; text-align: center; border-bottom-color: #FFFFFF; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3;" bgcolor="#FFFFFF" align="center">
      <td colspan="6" class="gt_heading gt_title gt_font_normal" style="border-style: none; color: #333333; font-size: 125%; padding-top: 4px; padding-bottom: 4px; padding-left: 5px; padding-right: 5px; border-bottom-width: 0; background-color: #FFFFFF; text-align: center; border-bottom-color: #FFFFFF; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; font-weight: normal;" bgcolor="#FFFFFF" align="center">S&amp;P 500</td>
    </tr>
    <tr class="gt_heading" style="border-style: none; background-color: #FFFFFF; text-align: center; border-bottom-color: #FFFFFF; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3;" bgcolor="#FFFFFF" align="center">
      <td colspan="6" class="gt_heading gt_subtitle gt_font_normal gt_bottom_border" style="border-style: none; color: #333333; font-size: 85%; padding-top: 3px; padding-bottom: 5px; padding-left: 5px; padding-right: 5px; border-top-color: #FFFFFF; border-top-width: 0; background-color: #FFFFFF; text-align: center; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; border-bottom-style: solid; border-bottom-width: 2px; border-bottom-color: #D3D3D3; font-weight: normal;" bgcolor="#FFFFFF" align="center">2010-06-07 to 2010-06-14</td>
    </tr>
    <tr class="gt_col_headings" style="border-style: none; border-top-style: solid; border-top-width: 2px; border-top-color: #D3D3D3; border-bottom-style: solid; border-bottom-width: 2px; border-bottom-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3;">
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="date" style="border-style: none; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: normal; text-transform: inherit; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: bottom; padding-top: 5px; padding-bottom: 6px; padding-left: 5px; padding-right: 5px; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" bgcolor="#FFFFFF" valign="bottom" align="right">date</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="open" style="border-style: none; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: normal; text-transform: inherit; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: bottom; padding-top: 5px; padding-bottom: 6px; padding-left: 5px; padding-right: 5px; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" bgcolor="#FFFFFF" valign="bottom" align="right">open</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="high" style="border-style: none; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: normal; text-transform: inherit; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: bottom; padding-top: 5px; padding-bottom: 6px; padding-left: 5px; padding-right: 5px; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" bgcolor="#FFFFFF" valign="bottom" align="right">high</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="low" style="border-style: none; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: normal; text-transform: inherit; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: bottom; padding-top: 5px; padding-bottom: 6px; padding-left: 5px; padding-right: 5px; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" bgcolor="#FFFFFF" valign="bottom" align="right">low</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="close" style="border-style: none; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: normal; text-transform: inherit; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: bottom; padding-top: 5px; padding-bottom: 6px; padding-left: 5px; padding-right: 5px; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" bgcolor="#FFFFFF" valign="bottom" align="right">close</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="volume" style="border-style: none; color: #333333; background-color: #FFFFFF; font-size: 100%; font-weight: normal; text-transform: inherit; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: bottom; padding-top: 5px; padding-bottom: 6px; padding-left: 5px; padding-right: 5px; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" bgcolor="#FFFFFF" valign="bottom" align="right">volume</th>
    </tr>
  </thead>
  <tbody class="gt_table_body" style="border-style: none; border-top-style: solid; border-top-width: 2px; border-top-color: #D3D3D3; border-bottom-style: solid; border-bottom-width: 2px; border-bottom-color: #D3D3D3;">
    <tr style="border-style: none;"><td headers="date" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">Mon, Jun 14, 2010</td>
<td headers="open" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">$1,095.00</td>
<td headers="high" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">$1,105.91</td>
<td headers="low" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">$1,089.03</td>
<td headers="close" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">$1,089.63</td>
<td headers="volume" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">4.43B</td></tr>
    <tr style="border-style: none;"><td headers="date" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">Fri, Jun 11, 2010</td>
<td headers="open" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">$1,082.65</td>
<td headers="high" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">$1,092.25</td>
<td headers="low" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">$1,077.12</td>
<td headers="close" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">$1,091.60</td>
<td headers="volume" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">4.06B</td></tr>
    <tr style="border-style: none;"><td headers="date" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">Thu, Jun 10, 2010</td>
<td headers="open" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">$1,058.77</td>
<td headers="high" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">$1,087.85</td>
<td headers="low" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">$1,058.77</td>
<td headers="close" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">$1,086.84</td>
<td headers="volume" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">5.14B</td></tr>
    <tr style="border-style: none;"><td headers="date" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">Wed, Jun 9, 2010</td>
<td headers="open" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">$1,062.75</td>
<td headers="high" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">$1,077.74</td>
<td headers="low" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">$1,052.25</td>
<td headers="close" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">$1,055.69</td>
<td headers="volume" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">5.98B</td></tr>
    <tr style="border-style: none;"><td headers="date" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">Tue, Jun 8, 2010</td>
<td headers="open" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">$1,050.81</td>
<td headers="high" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">$1,063.15</td>
<td headers="low" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">$1,042.17</td>
<td headers="close" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">$1,062.00</td>
<td headers="volume" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">6.19B</td></tr>
    <tr style="border-style: none;"><td headers="date" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">Mon, Jun 7, 2010</td>
<td headers="open" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">$1,065.84</td>
<td headers="high" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">$1,071.36</td>
<td headers="low" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">$1,049.86</td>
<td headers="close" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">$1,050.47</td>
<td headers="volume" class="gt_row gt_right" style="border-style: none; padding-top: 8px; padding-bottom: 8px; padding-left: 5px; padding-right: 5px; margin: 10px; border-top-style: solid; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: right; font-variant-numeric: tabular-nums;" valign="middle" align="right">5.47B</td></tr>
  </tbody>
  &#10;  
</table>
</div>

Show sequence alignment as a table - play around with this as an
(alternative to ggmsa):

This is a slice of SARS-CoV2’s Orf9b prot sequence aligned to related
viruses’ Orf9b

``` r
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

Show as a gt table, without using the as_raw_html() function

``` r
seq_slice_gt <- seq_slice_tbl %>% 
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

seq_slice_gt
```

<div id="ubjzvhcmuv" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>@import url("https://fonts.googleapis.com/css2?family=Courier:ital,wght@0,100;0,200;0,300;0,400;0,500;0,600;0,700;0,800;0,900;1,100;1,200;1,300;1,400;1,500;1,600;1,700;1,800;1,900&display=swap");
#ubjzvhcmuv table {
  font-family: Courier, system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#ubjzvhcmuv thead, #ubjzvhcmuv tbody, #ubjzvhcmuv tfoot, #ubjzvhcmuv tr, #ubjzvhcmuv td, #ubjzvhcmuv th {
  border-style: none;
}
&#10;#ubjzvhcmuv p {
  margin: 0;
  padding: 0;
}
&#10;#ubjzvhcmuv .gt_table {
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
&#10;#ubjzvhcmuv .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#ubjzvhcmuv .gt_title {
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
&#10;#ubjzvhcmuv .gt_subtitle {
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
&#10;#ubjzvhcmuv .gt_heading {
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
&#10;#ubjzvhcmuv .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#ubjzvhcmuv .gt_col_headings {
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
&#10;#ubjzvhcmuv .gt_col_heading {
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
&#10;#ubjzvhcmuv .gt_column_spanner_outer {
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
&#10;#ubjzvhcmuv .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#ubjzvhcmuv .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#ubjzvhcmuv .gt_column_spanner {
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
&#10;#ubjzvhcmuv .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#ubjzvhcmuv .gt_group_heading {
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
&#10;#ubjzvhcmuv .gt_empty_group_heading {
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
&#10;#ubjzvhcmuv .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#ubjzvhcmuv .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#ubjzvhcmuv .gt_row {
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
&#10;#ubjzvhcmuv .gt_stub {
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
&#10;#ubjzvhcmuv .gt_stub_row_group {
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
&#10;#ubjzvhcmuv .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#ubjzvhcmuv .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#ubjzvhcmuv .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#ubjzvhcmuv .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#ubjzvhcmuv .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#ubjzvhcmuv .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#ubjzvhcmuv .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#ubjzvhcmuv .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#ubjzvhcmuv .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#ubjzvhcmuv .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#ubjzvhcmuv .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: hidden;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#ubjzvhcmuv .gt_footnotes {
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
&#10;#ubjzvhcmuv .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#ubjzvhcmuv .gt_sourcenotes {
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
&#10;#ubjzvhcmuv .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#ubjzvhcmuv .gt_left {
  text-align: left;
}
&#10;#ubjzvhcmuv .gt_center {
  text-align: center;
}
&#10;#ubjzvhcmuv .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#ubjzvhcmuv .gt_font_normal {
  font-weight: normal;
}
&#10;#ubjzvhcmuv .gt_font_bold {
  font-weight: bold;
}
&#10;#ubjzvhcmuv .gt_font_italic {
  font-style: italic;
}
&#10;#ubjzvhcmuv .gt_super {
  font-size: 65%;
}
&#10;#ubjzvhcmuv .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#ubjzvhcmuv .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#ubjzvhcmuv .gt_indent_1 {
  text-indent: 5px;
}
&#10;#ubjzvhcmuv .gt_indent_2 {
  text-indent: 10px;
}
&#10;#ubjzvhcmuv .gt_indent_3 {
  text-indent: 15px;
}
&#10;#ubjzvhcmuv .gt_indent_4 {
  text-indent: 20px;
}
&#10;#ubjzvhcmuv .gt_indent_5 {
  text-indent: 25px;
}
&#10;#ubjzvhcmuv .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}
&#10;#ubjzvhcmuv div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
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

Now add as_raw_html() function

``` r
seq_slice_gt |>
    ## if I don't include as_raw_html(), when I render and display as a github page I get a lot of extra stuff above the table
    as_raw_html()
```

<div id="ksthkohlas" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
  &#10;  <table class="gt_table" style="-webkit-font-smoothing: antialiased; -moz-osx-font-smoothing: grayscale; font-family: Courier, system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji'; display: table; border-collapse: collapse; line-height: normal; margin-left: auto; margin-right: auto; color: #333333; font-size: 16px; font-weight: normal; font-style: normal; background-color: #FFFFFF; width: auto; border-top-style: solid; border-top-width: 2px; border-top-color: #A8A8A8; border-right-style: none; border-right-width: 2px; border-right-color: #D3D3D3; border-bottom-style: solid; border-bottom-width: 2px; border-bottom-color: #A8A8A8; border-left-style: none; border-left-width: 2px; border-left-color: #D3D3D3; table-layout: fixed;" data-quarto-disable-processing="false" data-quarto-bootstrap="false" bgcolor="#FFFFFF">
  <colgroup>
    <col style="width:75px;">
    <col>
    <col>
    <col>
    <col>
  </colgroup>
  &#10;  <tbody class="gt_table_body" style="border-style: none; border-top-style: solid; border-top-width: 2px; border-top-color: #D3D3D3; border-bottom-style: hidden; border-bottom-width: 2px; border-bottom-color: #D3D3D3;">
    <tr style="border-style: none;"><td headers="id" class="gt_row gt_left" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">seq1</td>
<td headers="pos1" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">F</td>
<td headers="pos2" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">Q</td>
<td headers="pos3" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">L</td>
<td headers="pos4" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">T</td></tr>
    <tr style="border-style: none;"><td headers="id" class="gt_row gt_left" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">seq2</td>
<td headers="pos1" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">F</td>
<td headers="pos2" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">Q</td>
<td headers="pos3" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">L</td>
<td headers="pos4" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">I</td></tr>
    <tr style="border-style: none;"><td headers="id" class="gt_row gt_left" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">seq3</td>
<td headers="pos1" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">F</td>
<td headers="pos2" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">Q</td>
<td headers="pos3" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">L</td>
<td headers="pos4" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">T</td></tr>
    <tr style="border-style: none;"><td headers="id" class="gt_row gt_left" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">seq4</td>
<td headers="pos1" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">F</td>
<td headers="pos2" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">R</td>
<td headers="pos3" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">L</td>
<td headers="pos4" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">T</td></tr>
    <tr style="border-style: none;"><td headers="id" class="gt_row gt_left" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">seq5</td>
<td headers="pos1" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">F</td>
<td headers="pos2" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">R</td>
<td headers="pos3" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">L</td>
<td headers="pos4" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">I</td></tr>
    <tr style="border-style: none;"><td headers="id" class="gt_row gt_left" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">seq6</td>
<td headers="pos1" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">F</td>
<td headers="pos2" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">R</td>
<td headers="pos3" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">L</td>
<td headers="pos4" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">T</td></tr>
    <tr style="border-style: none;"><td headers="id" class="gt_row gt_left" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">seq7</td>
<td headers="pos1" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">F</td>
<td headers="pos2" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">R</td>
<td headers="pos3" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">L</td>
<td headers="pos4" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">T</td></tr>
    <tr style="border-style: none;"><td headers="id" class="gt_row gt_left" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">seq8</td>
<td headers="pos1" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">F</td>
<td headers="pos2" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">Q</td>
<td headers="pos3" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">L</td>
<td headers="pos4" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">T</td></tr>
    <tr style="border-style: none;"><td headers="id" class="gt_row gt_left" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">seq9</td>
<td headers="pos1" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">F</td>
<td headers="pos2" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">Q</td>
<td headers="pos3" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">S</td>
<td headers="pos4" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">T</td></tr>
    <tr style="border-style: none;"><td headers="id" class="gt_row gt_left" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: left;" valign="middle" align="left">seq10</td>
<td headers="pos1" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">F</td>
<td headers="pos2" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">Q</td>
<td headers="pos3" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">S</td>
<td headers="pos4" class="gt_row gt_center" style="border-style: none; padding-top: 0px; padding-bottom: 0px; padding-left: 2px; padding-right: 2px; margin: 10px; border-top-style: hidden; border-top-width: 1px; border-top-color: #D3D3D3; border-left-style: none; border-left-width: 1px; border-left-color: #D3D3D3; border-right-style: none; border-right-width: 1px; border-right-color: #D3D3D3; vertical-align: middle; overflow-x: hidden; text-align: center;" valign="middle" align="center">T</td></tr>
  </tbody>
  &#10;  
</table>
</div>

# Finished

``` r
sessionInfo()
```

    ## R version 4.5.1 (2025-06-13)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Sequoia 15.6.1
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
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
    ##  [5] stringr_1.5.1   dplyr_1.1.4     purrr_1.1.0     readr_2.1.5    
    ##  [9] tidyr_1.3.1     tibble_3.3.0    ggplot2_3.5.2   tidyverse_2.0.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] jsonlite_2.0.0     gtable_0.3.6       compiler_4.5.1     Rcpp_1.1.0        
    ##  [5] tidyselect_1.2.1   xml2_1.3.8         scales_1.4.0       yaml_2.3.10       
    ##  [9] fastmap_1.2.0      R6_2.6.1           generics_0.1.4     curl_6.4.0        
    ## [13] bigD_0.3.1         knitr_1.50         rprojroot_2.1.0    pillar_1.11.0     
    ## [17] RColorBrewer_1.1-3 tzdb_0.5.0         rlang_1.1.6        V8_7.0.0          
    ## [21] stringi_1.8.7      xfun_0.52          sass_0.4.10        juicyjuice_0.1.0  
    ## [25] timechange_0.3.0   cli_3.6.5          withr_3.0.2        magrittr_2.0.3    
    ## [29] digest_0.6.37      grid_4.5.1         rstudioapi_0.17.1  hms_1.1.3         
    ## [33] lifecycle_1.0.4    vctrs_0.6.5        evaluate_1.0.4     glue_1.8.0        
    ## [37] farver_2.1.2       rmarkdown_2.29     tools_4.5.1        pkgconfig_2.0.3   
    ## [41] htmltools_0.5.8.1

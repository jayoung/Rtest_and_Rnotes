table_display_patchwork_and_gt_demo
================
Janet Young

2025-11-07

Testing a few details to do with wrapping long column names when I
display tables, as well as combining tables and plots using the
`patchwork` package.

Make an example tibble with long names (`my_tbl`):

``` r
my_tbl <- tibble(
    first_column_break_with_a_very_break_long_name = c("apple", "orange"),
    second_column_break_with_a_long_name = c(1034,1365)
)
my_tbl
```

    ## # A tibble: 2 × 2
    ##   first_column_break_with_a_very_break_long_name second_column_break_with_a_lo…¹
    ##   <chr>                                                                    <dbl>
    ## 1 apple                                                                     1034
    ## 2 orange                                                                    1365
    ## # ℹ abbreviated name: ¹​second_column_break_with_a_long_name

Create `my_tbl_gt` - I use gt to display that, inserting `<br>`
linebreaks, and doing some formatting (italics, and font size). Looks
good.

``` r
my_tbl_gt <- gt(my_tbl) %>% 
    fmt_number(use_seps = TRUE, decimals=0) %>% 
    cols_label(first_column_break_with_a_very_break_long_name="<em>first column</em><br>with a very<br>long name",
               second_column_break_with_a_long_name="second column<br>with a long name",
               .fn = md ) %>%  
    ## this makes font 75% the size of default
    tab_style(style = cell_text(size = pct(75)),
              locations=list(cells_column_labels(), cells_body()))
my_tbl_gt
```

<div id="bqsfreeusd" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#bqsfreeusd table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#bqsfreeusd thead, #bqsfreeusd tbody, #bqsfreeusd tfoot, #bqsfreeusd tr, #bqsfreeusd td, #bqsfreeusd th {
  border-style: none;
}
&#10;#bqsfreeusd p {
  margin: 0;
  padding: 0;
}
&#10;#bqsfreeusd .gt_table {
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
&#10;#bqsfreeusd .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#bqsfreeusd .gt_title {
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
&#10;#bqsfreeusd .gt_subtitle {
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
&#10;#bqsfreeusd .gt_heading {
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
&#10;#bqsfreeusd .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#bqsfreeusd .gt_col_headings {
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
&#10;#bqsfreeusd .gt_col_heading {
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
&#10;#bqsfreeusd .gt_column_spanner_outer {
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
&#10;#bqsfreeusd .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#bqsfreeusd .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#bqsfreeusd .gt_column_spanner {
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
&#10;#bqsfreeusd .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#bqsfreeusd .gt_group_heading {
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
&#10;#bqsfreeusd .gt_empty_group_heading {
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
&#10;#bqsfreeusd .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#bqsfreeusd .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#bqsfreeusd .gt_row {
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
&#10;#bqsfreeusd .gt_stub {
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
&#10;#bqsfreeusd .gt_stub_row_group {
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
&#10;#bqsfreeusd .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#bqsfreeusd .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#bqsfreeusd .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#bqsfreeusd .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#bqsfreeusd .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#bqsfreeusd .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#bqsfreeusd .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#bqsfreeusd .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#bqsfreeusd .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#bqsfreeusd .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#bqsfreeusd .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#bqsfreeusd .gt_footnotes {
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
&#10;#bqsfreeusd .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#bqsfreeusd .gt_sourcenotes {
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
&#10;#bqsfreeusd .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#bqsfreeusd .gt_left {
  text-align: left;
}
&#10;#bqsfreeusd .gt_center {
  text-align: center;
}
&#10;#bqsfreeusd .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#bqsfreeusd .gt_font_normal {
  font-weight: normal;
}
&#10;#bqsfreeusd .gt_font_bold {
  font-weight: bold;
}
&#10;#bqsfreeusd .gt_font_italic {
  font-style: italic;
}
&#10;#bqsfreeusd .gt_super {
  font-size: 65%;
}
&#10;#bqsfreeusd .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#bqsfreeusd .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#bqsfreeusd .gt_indent_1 {
  text-indent: 5px;
}
&#10;#bqsfreeusd .gt_indent_2 {
  text-indent: 10px;
}
&#10;#bqsfreeusd .gt_indent_3 {
  text-indent: 15px;
}
&#10;#bqsfreeusd .gt_indent_4 {
  text-indent: 20px;
}
&#10;#bqsfreeusd .gt_indent_5 {
  text-indent: 25px;
}
&#10;#bqsfreeusd .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}
&#10;#bqsfreeusd div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" style="font-size: 75%;" scope="col" id="first_column_break_with_a_very_break_long_name"><span class='gt_from_md'><em>first column</em><br>with a very<br>long name</span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" style="font-size: 75%;" scope="col" id="second_column_break_with_a_long_name"><span class='gt_from_md'>second column<br>with a long name</span></th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="first_column_break_with_a_very_break_long_name" class="gt_row gt_left" style="font-size: 75%;">apple</td>
<td headers="second_column_break_with_a_long_name" class="gt_row gt_right" style="font-size: 75%;">1,034</td></tr>
    <tr><td headers="first_column_break_with_a_very_break_long_name" class="gt_row gt_left" style="font-size: 75%;">orange</td>
<td headers="second_column_break_with_a_long_name" class="gt_row gt_right" style="font-size: 75%;">1,365</td></tr>
  </tbody>
  &#10;  
</table>
</div>

Now display `my_tbl_gt` using `patchwork::wrap_table()`. It still wraps
the colnames, and keeps the column widths reasonable (although it does
drop the italics):

``` r
options(gt.html_tag_check = FALSE)
wrap_table(my_tbl_gt, panel="full")
```

![](table_display_patchwork_and_gt_demo_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

We can also give the table a caption using gt. Looks OK if I display
using gt:

``` r
my_tbl_gt_with_long_caption <- gt(my_tbl, caption="My caption which is very very long with a lot of words") %>% 
    fmt_number(use_seps = TRUE, decimals=0) %>% 
    cols_label(first_column_break_with_a_very_break_long_name="<em>first column</em><br>with a very<br>long name",
               second_column_break_with_a_long_name="second column<br>with a long name",
               .fn = md ) %>%  
    ## this makes font 75% the size of default
    tab_style(style = cell_text(size = pct(75)),
              locations=list(cells_column_labels(),
                             cells_body()))
my_tbl_gt_with_long_caption
```

<div id="egssyygwyj" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#egssyygwyj table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#egssyygwyj thead, #egssyygwyj tbody, #egssyygwyj tfoot, #egssyygwyj tr, #egssyygwyj td, #egssyygwyj th {
  border-style: none;
}
&#10;#egssyygwyj p {
  margin: 0;
  padding: 0;
}
&#10;#egssyygwyj .gt_table {
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
&#10;#egssyygwyj .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#egssyygwyj .gt_title {
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
&#10;#egssyygwyj .gt_subtitle {
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
&#10;#egssyygwyj .gt_heading {
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
&#10;#egssyygwyj .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#egssyygwyj .gt_col_headings {
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
&#10;#egssyygwyj .gt_col_heading {
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
&#10;#egssyygwyj .gt_column_spanner_outer {
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
&#10;#egssyygwyj .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#egssyygwyj .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#egssyygwyj .gt_column_spanner {
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
&#10;#egssyygwyj .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#egssyygwyj .gt_group_heading {
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
&#10;#egssyygwyj .gt_empty_group_heading {
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
&#10;#egssyygwyj .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#egssyygwyj .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#egssyygwyj .gt_row {
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
&#10;#egssyygwyj .gt_stub {
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
&#10;#egssyygwyj .gt_stub_row_group {
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
&#10;#egssyygwyj .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#egssyygwyj .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#egssyygwyj .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#egssyygwyj .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#egssyygwyj .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#egssyygwyj .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#egssyygwyj .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#egssyygwyj .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#egssyygwyj .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#egssyygwyj .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#egssyygwyj .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#egssyygwyj .gt_footnotes {
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
&#10;#egssyygwyj .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#egssyygwyj .gt_sourcenotes {
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
&#10;#egssyygwyj .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#egssyygwyj .gt_left {
  text-align: left;
}
&#10;#egssyygwyj .gt_center {
  text-align: center;
}
&#10;#egssyygwyj .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#egssyygwyj .gt_font_normal {
  font-weight: normal;
}
&#10;#egssyygwyj .gt_font_bold {
  font-weight: bold;
}
&#10;#egssyygwyj .gt_font_italic {
  font-style: italic;
}
&#10;#egssyygwyj .gt_super {
  font-size: 65%;
}
&#10;#egssyygwyj .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#egssyygwyj .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#egssyygwyj .gt_indent_1 {
  text-indent: 5px;
}
&#10;#egssyygwyj .gt_indent_2 {
  text-indent: 10px;
}
&#10;#egssyygwyj .gt_indent_3 {
  text-indent: 15px;
}
&#10;#egssyygwyj .gt_indent_4 {
  text-indent: 20px;
}
&#10;#egssyygwyj .gt_indent_5 {
  text-indent: 25px;
}
&#10;#egssyygwyj .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}
&#10;#egssyygwyj div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <caption>My caption which is very very long with a lot of words</caption>
  <thead>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" style="font-size: 75%;" scope="col" id="first_column_break_with_a_very_break_long_name"><span class='gt_from_md'><em>first column</em><br>with a very<br>long name</span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" style="font-size: 75%;" scope="col" id="second_column_break_with_a_long_name"><span class='gt_from_md'>second column<br>with a long name</span></th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="first_column_break_with_a_very_break_long_name" class="gt_row gt_left" style="font-size: 75%;">apple</td>
<td headers="second_column_break_with_a_long_name" class="gt_row gt_right" style="font-size: 75%;">1,034</td></tr>
    <tr><td headers="first_column_break_with_a_very_break_long_name" class="gt_row gt_left" style="font-size: 75%;">orange</td>
<td headers="second_column_break_with_a_long_name" class="gt_row gt_right" style="font-size: 75%;">1,365</td></tr>
  </tbody>
  &#10;  
</table>
</div>

However, if we display a gt table with a long caption using
`patchwork::wrap_table()` (which we might do if we want to combine it
with plots), then two things are wrong:

- the caption is not visible (a knowon issue - we could deal with that
  by adding `+ plot_annotation(title="my title")`).
- worse: even though the caption is not visible, its presence messes up
  the column widths:

``` r
wrap_table(my_tbl_gt_with_long_caption, panel="full")
```

![](table_display_patchwork_and_gt_demo_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

We can also make a version where we include `\n` in the colnames, to try
to wrap them a different way

``` r
my_tbl_wrap_colnames <- my_tbl %>% 
    set_names(nm=c("first column\nwith a very\nlong name",
                   "second column\nwith a long name"))
```

A plain print output does not wrap colnames:

``` r
my_tbl_wrap_colnames
```

    ## # A tibble: 2 × 2
    ##   `first column\nwith a very\nlong name` `second column\nwith a long name`
    ##   <chr>                                                              <dbl>
    ## 1 apple                                                               1034
    ## 2 orange                                                              1365

Neither does using `gt()`:

``` r
my_tbl_wrap_colnames %>% gt()
```

<div id="nlqgitariz" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#nlqgitariz table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#nlqgitariz thead, #nlqgitariz tbody, #nlqgitariz tfoot, #nlqgitariz tr, #nlqgitariz td, #nlqgitariz th {
  border-style: none;
}
&#10;#nlqgitariz p {
  margin: 0;
  padding: 0;
}
&#10;#nlqgitariz .gt_table {
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
&#10;#nlqgitariz .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#nlqgitariz .gt_title {
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
&#10;#nlqgitariz .gt_subtitle {
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
&#10;#nlqgitariz .gt_heading {
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
&#10;#nlqgitariz .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#nlqgitariz .gt_col_headings {
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
&#10;#nlqgitariz .gt_col_heading {
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
&#10;#nlqgitariz .gt_column_spanner_outer {
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
&#10;#nlqgitariz .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#nlqgitariz .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#nlqgitariz .gt_column_spanner {
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
&#10;#nlqgitariz .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#nlqgitariz .gt_group_heading {
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
&#10;#nlqgitariz .gt_empty_group_heading {
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
&#10;#nlqgitariz .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#nlqgitariz .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#nlqgitariz .gt_row {
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
&#10;#nlqgitariz .gt_stub {
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
&#10;#nlqgitariz .gt_stub_row_group {
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
&#10;#nlqgitariz .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#nlqgitariz .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#nlqgitariz .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#nlqgitariz .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#nlqgitariz .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#nlqgitariz .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#nlqgitariz .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#nlqgitariz .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#nlqgitariz .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#nlqgitariz .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#nlqgitariz .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#nlqgitariz .gt_footnotes {
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
&#10;#nlqgitariz .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#nlqgitariz .gt_sourcenotes {
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
&#10;#nlqgitariz .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#nlqgitariz .gt_left {
  text-align: left;
}
&#10;#nlqgitariz .gt_center {
  text-align: center;
}
&#10;#nlqgitariz .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#nlqgitariz .gt_font_normal {
  font-weight: normal;
}
&#10;#nlqgitariz .gt_font_bold {
  font-weight: bold;
}
&#10;#nlqgitariz .gt_font_italic {
  font-style: italic;
}
&#10;#nlqgitariz .gt_super {
  font-size: 65%;
}
&#10;#nlqgitariz .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#nlqgitariz .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#nlqgitariz .gt_indent_1 {
  text-indent: 5px;
}
&#10;#nlqgitariz .gt_indent_2 {
  text-indent: 10px;
}
&#10;#nlqgitariz .gt_indent_3 {
  text-indent: 15px;
}
&#10;#nlqgitariz .gt_indent_4 {
  text-indent: 20px;
}
&#10;#nlqgitariz .gt_indent_5 {
  text-indent: 25px;
}
&#10;#nlqgitariz .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}
&#10;#nlqgitariz div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="first-column-with-a-very-long-name">first column
with a very
long name</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="second-column-with-a-long-name">second column
with a long name</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="first column
with a very
long name" class="gt_row gt_left">apple</td>
<td headers="second column
with a long name" class="gt_row gt_right">1034</td></tr>
    <tr><td headers="first column
with a very
long name" class="gt_row gt_left">orange</td>
<td headers="second column
with a long name" class="gt_row gt_right">1365</td></tr>
  </tbody>
  &#10;  
</table>
</div>

But `wrap_table()` does see the `\n` in colnames:

``` r
wrap_table(my_tbl_wrap_colnames)
```

![](table_display_patchwork_and_gt_demo_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

Combining the table with a plot (without gt):

(harder to add styling like font size this way)

``` r
p1 <- tibble(x=1:10, y=1:10) %>% 
    ggplot(aes(x=x,y=y))+
    geom_point() +
    theme_classic()

(p1 + wrap_table(my_tbl_wrap_colnames, panel="full") + p1) +
    plot_layout(widths=c(1,2,1))
```

![](table_display_patchwork_and_gt_demo_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

Combining the table with a plot (WITH gt):

``` r
p1 <- tibble(x=1:10, y=1:10) %>% 
    ggplot(aes(x=x,y=y))+
    geom_point() +
    theme_classic()

(p1 + wrap_table(my_tbl_gt, panel="full") + p1) +
    plot_layout(widths=c(1,2,1))
```

![](table_display_patchwork_and_gt_demo_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
sessionInfo()
```

    ## R version 4.5.1 (2025-06-13)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Sequoia 15.7.1
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
    ##  [1] gt_1.0.0        patchwork_1.3.2 lubridate_1.9.4 forcats_1.0.0  
    ##  [5] stringr_1.5.2   dplyr_1.1.4     purrr_1.1.0     readr_2.1.5    
    ##  [9] tidyr_1.3.1     tibble_3.3.0    ggplot2_3.5.2   tidyverse_2.0.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.6       compiler_4.5.1     tidyselect_1.2.1   xml2_1.4.0        
    ##  [5] scales_1.4.0       yaml_2.3.10        fastmap_1.2.0      R6_2.6.1          
    ##  [9] labeling_0.4.3     commonmark_2.0.0   generics_0.1.4     knitr_1.50        
    ## [13] pillar_1.11.1      RColorBrewer_1.1-3 tzdb_0.5.0         rlang_1.1.6       
    ## [17] utf8_1.2.6         stringi_1.8.7      litedown_0.7       xfun_0.53         
    ## [21] sass_0.4.10        timechange_0.3.0   cli_3.6.5          withr_3.0.2       
    ## [25] magrittr_2.0.4     digest_0.6.37      grid_4.5.1         rstudioapi_0.17.1 
    ## [29] markdown_2.0       hms_1.1.3          lifecycle_1.0.4    vctrs_0.6.5       
    ## [33] evaluate_1.0.5     glue_1.8.0         farver_2.1.2       rmarkdown_2.29    
    ## [37] tools_4.5.1        pkgconfig_2.0.3    htmltools_0.5.8.1

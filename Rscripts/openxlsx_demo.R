library(openxlsx)

# Tim's example (Hutch Slack channel)
wb <- createWorkbook()
addWorksheet(wb, "SheetName")
freezePane(wb, "SheetName", firstRow = TRUE)
writeDataTable(wb, "SheetName", iris, tableStyle = "TableStyleLight9")
saveWorkbook(wb, "Desktop/temp.xlsx", overwrite = TRUE)

# with my addition
myStyle <- createStyle(fontName="Arial", fontSize = 12, textDecoration="bold")
addStyle(wb, "SheetName", myStyle, stack=TRUE, 
         cols=c(2,3), rows=1:dim(iris)[1], gridExpand=TRUE)
saveWorkbook(wb, "Desktop/temp.xlsx", overwrite = TRUE)


# example from https://cran.r-project.org/web/packages/openxlsx/vignettes/formatting.pdf

library(openxlxs)



## data.frame to write
df <- data.frame("Date" = Sys.Date()-0:4,
                 "Logical" = c(TRUE, FALSE, TRUE, TRUE, FALSE),
                 "Currency" = paste("$",-2:2),
                 "Accounting" = -2:2,
                 "hLink" = "https://CRAN.R-project.org/",
                 "Percentage" = seq(-1, 1, length.out=5),
                 "TinyNumber" = runif(5) / 1E9, stringsAsFactors = FALSE)
class(df$Currency) <- "currency"
class(df$Accounting) <- "accounting"
class(df$hLink) <- "hyperlink"
class(df$Percentage) <- "percentage"
class(df$TinyNumber) <- "scientific"


## Formatting can be applied simply through the write functions
## global options can be set to further simplify things
options("openxlsx.borderStyle" = "thin")
options("openxlsx.borderColour" = "#4F81BD")

## create a workbook and add a worksheet
wb <- createWorkbook()
addWorksheet(wb, "writeData auto-formatting")
writeData(wb, 1, df, startRow = 2, startCol = 2) # df goes into the sheet starting at cell B2
saveWorkbook(wb, "Desktop/temp.xlsx", overwrite = TRUE)

writeData(wb, 1, df, startRow = 9, startCol = 2, borders = "surrounding") # df goes into the sheet again 
saveWorkbook(wb, "Desktop/temp.xlsx", overwrite = TRUE)

writeData(wb, 1, df, startRow = 16, startCol = 2, borders = "rows")
writeData(wb, 1, df, startRow = 23, startCol = 2, borders ="columns")
writeData(wb, 1, df, startRow = 30, startCol = 2, borders ="all")
saveWorkbook(wb, "Desktop/temp.xlsx", overwrite = TRUE)


## headerStyles
hs1 <- createStyle(fgFill = "#4F81BD", halign = "CENTER", 
                   textDecoration = "Bold",
                   border = "Bottom", fontColour = "white")
writeData(wb, 1, df, startRow = 16, startCol = 10, 
          headerStyle = hs1,
          borders = "rows", borderStyle = "medium")
saveWorkbook(wb, "Desktop/temp.xlsx", overwrite = TRUE)

## to change the display text for a hyperlink column just write over those cells
writeData(wb, sheet = 1, x = paste("Hyperlink", 1:5), startRow = 17, startCol = 14)
saveWorkbook(wb, "Desktop/temp.xlsx", overwrite = TRUE)

## writing as an Excel Table
addWorksheet(wb, "writeDataTable")
writeDataTable(wb, 2, df, startRow = 2, startCol = 2)
# TableStyleLight9 etc are defined within Excel
writeDataTable(wb, 2, df, startRow = 9, startCol = 2, tableStyle = "TableStyleLight9")
writeDataTable(wb, 2, df, startRow = 16, startCol = 2, tableStyle = "TableStyleLight2")
writeDataTable(wb, 2, df, startRow = 23, startCol = 2, tableStyle = "TableStyleMedium21")
openXL(wb) ## opens a temp version - pops up an Excel window

saveWorkbook(wb, "Desktop/temp2.xlsx", overwrite = TRUE)



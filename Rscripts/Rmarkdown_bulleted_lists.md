Rmarkdown bulleted lists
================
Janet Young

2025-02-10

Rules for bulleted lists:

- there should be an EMPTY LINE before the list starts
- top-level list items start with an `*`, `-`, or `+` and don’t need an
  indent. I think I’ll use `*` to keep things uniform.
  - next-level list items are indented. Can also start with `*`, `-`, or
    `+` but I think I’ll use `+` to keep things uniform.
    - third-level indent
- you can add more items

Rules for numbered lists:

1.  there should be an EMPTY LINE before the list starts
2.  top-level list items start with an asterisk and don’t need an indent
    - next-level list items start with a letter
    - second next-level item here
    - it doesn’t seem easy to have a/b/c tags for second level. Maybe
      [this link](https://pandoc.org/MANUAL.html#ordered-lists) explains
      how, but it’s unclear.
3.  you can add more items

Without the blank line before the list, you might be able to make it
work somehow, but you need lots of extra spaces at the end of pretty
much every line

Here’s a bulleted list WITHOUT the empty line and WITHOUT extra spaces -
it looks really bad: \* there should be an EMPTY LINE before the list
starts \* top-level list items start with an `*`, `-`, or `+` and don’t
need an indent. I think I’ll use `*` to keep things uniform. +
next-level list items are indented. Can also start with `*`, `-`, or `+`
but I think I’ll use `+` to keep things uniform. - third-level indent \*
you can add more items

Here’s a bulleted list WITHOUT the empty line but WITH extra spaces. It
still looks bad!  
\* there should be an EMPTY LINE before the list starts.  
\* top-level list items start with an `*`, `-`, or `+` and don’t need an
indent. I think I’ll use `*` to keep things uniform.  
+ next-level list items are indented. Can also start with `*`, `-`, or
`+` but I think I’ll use `+` to keep things uniform.  
- third-level indent.  
\* you can add more items.

library(rdrop2)
# Authenticate and save token for later use2
token <- drop_auth(rdstoken = "dropbox_token.rds")

# Retrieveing your file is as simple as
drop_download("public/conclusions/IGHV1-18.docx", local_path = "IGHV1-18.docx",
              overwrite = TRUE)
drop_df <- textreadr::read_docx("IGHV1-18.docx")

system("pandoc -s IGHV1-18.docx --wrap=none --reference-links -t markdown -o IGHV1-18_dock.Rmd")

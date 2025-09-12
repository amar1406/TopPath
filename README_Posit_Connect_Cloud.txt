Posit Connect Cloud (Free) — Quick Deploy Steps

1) Put these files into your project at `ShinnyApp/`:
   - app.R  (this wrapper)
   - requirements.txt
   - .Rprofile

2) Remove Windows-specific paths from your scripts:
   - Delete any setwd("C:/...") calls.
   - Use file.path("database", "your_file") for data access.

3) Initialize renv locally (R console, inside ShinnyApp/):
   install.packages("renv")
   renv::init()
   renv::snapshot()

4) Deploy to Posit Connect Cloud from R:
   install.packages("rsconnect")
   rsconnect::deployApp("ShinnyApp")

5) After deploy, open the app and verify:
   - It loads without errors.
   - Server logs should show py_config() with Python + torch versions.

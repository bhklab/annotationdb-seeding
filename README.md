# AnnotationDB Seeding

Repository for building and populating tables in the AnnotationDB MySQL database. Specifically, this database coordinates the creation and population of a drug database for PubChem, a ChEMBL drug database with FKs to pubchem CIDs, and a cell lines database using cellosaurus. This repository utilizes Pixi for package management and SQLAlchemy for main database operations.

## To Run Scripts


Environment Variables
```Bash
DATABASE_IP=
DATABASE_PASS=
DATABASE_USER=
```
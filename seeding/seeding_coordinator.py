from sqlalchemy import create_engine
import os
from dotenv import load_dotenv
from urllib.parse import quote_plus
import pandas as pd

from create_tables import Base, Compounds, CompoundSynonyms, ChemblMechanism, CellLines

load_dotenv(override=True)

password_cleaned = quote_plus(os.getenv("DATABASE_PASS"))
engine = create_engine(
    f"mysql+pymysql://{os.getenv('DATABASE_USER')}:{password_cleaned}"
    f"@{os.getenv('DATABASE_IP')}:{os.getenv('PORT')}/{os.getenv('SELECTED_DB')}",
    echo=True,
)

Base.metadata.create_all(engine)


# Helper to align a dataframe to a model's columns
def align_to_model(df: pd.DataFrame, model) -> pd.DataFrame:
    cols = [c.name for c in model.__table__.columns]
    return df[[c for c in cols if c in df.columns]]


# Read the attached files directly
compounds_df = pd.read_csv("data_extraction/drugs/pubchem/output_data/drug_out.csv")
synonyms_df = pd.read_csv("data_extraction/drugs/pubchem/output_data/drug_synonyms.csv")
chembl_mech_df = pd.read_csv(
    "data_extraction/drugs/pubchem/output_data/chembl_mechanism.csv"
)
cell_lines_df = pd.read_csv(
    "data_extraction/cell_lines/cellosaurus/output_data/cell_lines_table_cleaned.csv"
)

# Align columns to ORM models
compounds_df = align_to_model(compounds_df, Compounds)
synonyms_df = align_to_model(synonyms_df, CompoundSynonyms)
chembl_mech_df = align_to_model(chembl_mech_df, ChemblMechanism)
cell_lines_df = align_to_model(cell_lines_df, CellLines)

# Insert into tables named by the ORM models
compounds_df.to_sql(
    name=Compounds.__tablename__, con=engine, if_exists="append", index=False
)
synonyms_df.to_sql(
    name=CompoundSynonyms.__tablename__, con=engine, if_exists="append", index=False
)
chembl_mech_df.to_sql(
    name=ChemblMechanism.__tablename__, con=engine, if_exists="append", index=False
)
cell_lines_df.to_sql(
    name=CellLines.__tablename__, con=engine, if_exists="append", index=False
)

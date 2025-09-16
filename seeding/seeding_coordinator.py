from sqlalchemy import create_engine
import os
from dotenv import load_dotenv
from create_tables import Base
from urllib.parse import quote_plus
import pandas as pd


load_dotenv(override=True)

password_cleaned = quote_plus(os.getenv("DATABASE_PASS"))
engine = create_engine(
    f"mysql+pymysql://{os.getenv('DATABASE_USER')}:{password_cleaned}"
    f"@{os.getenv('DATABASE_IP')}:{os.getenv('PORT')}/{os.getenv('SELECTED_DB')}",
    echo=True,
)

Base.metadata.create_all(engine)


df = pd.read_csv(
    "../data_extraction/drugs/pubchem/output_data/pubchem_table_omitted.csv"
)
df = df.drop(columns=[df.columns[0]])  # Remove first empty column

df.to_sql(name="pubchem", con=engine, if_exists="append", index=False)

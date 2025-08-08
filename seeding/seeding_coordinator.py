import os
from urllib.parse import quote_plus
from sqlalchemy import create_engine
from sqlalchemy.orm import Session
from create_tables import Base
from dotenv import load_dotenv

load_dotenv(override=True)

DATABASE_URL=f"mysql+pymysql://{os.getenv("DATABASE_USER")}:{quote_plus(os.getenv("DATABASE_PASS"))}@{os.getenv("DATABASE_IP")}:{os.getenv("PORT")}/{os.getenv("SELECTED_DB")}"

# Database engine (connection setup)
engine = create_engine(DATABASE_URL, echo=True)

# Create tables
Base.metadata.create_all(engine)

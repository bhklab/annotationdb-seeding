from sqlalchemy import create_engine
from sqlalchemy.orm import Session
from create_tables import Pubchem, Chembl, Cellosaurus
import os
from dotenv import load_dotenv

load_dotenv(override=True)

engine = create_engine(os.getenv("DATABASE_IP"), echo=True)

Base.metadata.create_all(engine)

with Session(engine) as Session:

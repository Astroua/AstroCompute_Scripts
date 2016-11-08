from app import db
from app.models import Data

db.create_all()

print("DB created.")
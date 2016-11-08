from app import db

class Data(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    nickname = db.Column(db.String(64), index=True, unique=True)
    password = db.Column(db.String(120), index=True, unique=True)
    inst = db.Column(db.String(120), index=True, unique=True)
    emailaddress = db.Column(db.String(120), index=True, unique=True)

    @classmethod
    def query_by_name(cls, user_name,pass_word):
	user_query=Data.query.filter_by(nickname=user_name.encode('ascii','ignore')).first()
	if user_query!=None:
		if user_query.password.encode('ascii','ignore')==pass_word:
			return([user_query.nickname,user_query.password])
		else:
			return([user_query.nickname,'no_password'])
	else:
		return('no_user')

    def __repr__(self):
        return '<Data %r>' % (self.nickname)
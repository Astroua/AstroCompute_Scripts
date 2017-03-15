from app import db

from user import User


class Data(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    nickname = db.Column(db.String(64), index=True, unique=True)
    password = db.Column(db.String(120), index=True, unique=True)
    inst = db.Column(db.String(120), index=True, unique=True)
    emailaddress = db.Column(db.String(120), index=True, unique=True)

    @classmethod
    def query_by_name(cls, user_name, pass_word):
        user_query = Data.query.filter_by(
            nickname=user_name.encode('ascii', 'ignore')).first()
        if user_query is not None:
            if user_query.password.encode('ascii', 'ignore') == pass_word:
                return(User(user_query))
            else:
                return('no_password')
        else:
            return('no_user')

    @classmethod
    def user_query(cls, user_name):
        '''
        Return user query if the user_name is matched in the database.
        '''
        user_query = Data.query.filter_by(
            nickname=user_name.encode('ascii', 'ignore')).first()
        if user_query is not None:
                return(User(user_query))

    def __repr__(self):
        return '<Data %r>' % (self.nickname)

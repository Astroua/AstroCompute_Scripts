
from flask.ext.login import UserMixin


class User(UserMixin):

    def __init__(self, query_result, active=True):
        self.username = query_result.nickname
        self.password = query_result.password
        self.institution = query_result.inst
        self.email = query_result.emailaddress

        self.id = query_result.nickname

        self.active = active

    def is_active(self):
        # Here you should write whatever the code is
        # that checks the database if your user is active
        return self.active

    def is_authenticated(self):
        return True

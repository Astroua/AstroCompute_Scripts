
from wtforms import StringField
from wtforms.widgets import html_params, HTMLString


class FolderInput(object):
    """
    Renders a file input chooser field.
    """

    def __call__(self, field, **kwargs):
        kwargs.setdefault('id', field.id)
        return HTMLString('<input %s id="ctrl" webkitdirectory directory multiple>' % \
                          html_params(name=field.name, type='file', **kwargs))


class FolderField(StringField):
    """
    Can render a file-upload field.  Will take any passed filename value, if
    any is sent by the browser in the post params.  This field will NOT
    actually handle the file upload portion, as wtforms does not deal with
    individual frameworks' file handling capabilities.
    """
    widget = FolderInput()

import os
import jinja2
import subprocess


def makeDirIfNotExists(filePath):
    '''if no directory for files create it'''
    if not os.path.exists(os.path.dirname(filePath)):
        try:
            os.makedirs(os.path.dirname(filePath))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise


class Templater(object):

    def __init__( self, outputFilePath, templatePath, values ):

        self.outputFilePath = outputFilePath
        self.templatePath = templatePath
        self.values = values
        
        self.create_text_file(self.outputFilePath, self.templatePath, self.values)
        

    def template_engine(self, template, values):

        '''creates a text from a template'''
    
        # Jinja environment must be set to be compatible with other syntaxes
        jinja_env = jinja2.Environment(
	          block_start_string = '{%',
	          block_end_string = '%}',
	          variable_start_string = '{{',
	          variable_end_string = '}}',
	          comment_start_string = '!#',
	          comment_end_string = '#!',
	          line_statement_prefix = '%%',
	          line_comment_prefix = '%#',
	          trim_blocks = True,
	          autoescape = False,
	          loader = jinja2.FileSystemLoader('./')
        )

        # define which template to use
        template = jinja_env.get_template(template) 
        
        # input values into a template and create text
        text = template.render(values)

        return text

    
    def create_text_file(self, file_path, template, values):
        '''creates a file from a template'''
        
        makeDirIfNotExists(self.outputFilePath)
        
        # create text
        text = self.template_engine(template,values)
        
        # create a file
        out_file = open(file_path, 'w+')
        out_file.write(text)
        out_file.close()
        
        return file_path

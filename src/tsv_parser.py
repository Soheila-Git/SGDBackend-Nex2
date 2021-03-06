import csv
import os
import json
from sqlalchemy.exc import IntegrityError
import traceback

from loading.load_summaries_sync import load_summaries
from helpers import upload_file

# takes a TSV file and returns an array of annotations
def parse_tsv_annotations(db_session, tsv_file, filename, template_type, username):
    db_session.execute('SET LOCAL ROLE ' + username)
    try:
        if not filename.endswith('.tsv'):
            raise ValueError('File format not accepted. Please upload a valid TSV file.')
        raw_file_content = csv.reader(tsv_file, delimiter='\t', dialect=csv.excel_tab)
        tsv_file.seek(0)
    except:
        traceback.print_exc()
        db_session.close()
        raise ValueError('File format not accepted. Please upload a valid TSV file.')
    try:
	    upload_file(
            username, tsv_file,
            filename=filename,
            data_id=248375,
            description='summary upload',
            display_name=filename,
            format_id=248824,
            format_name='TSV',
            file_extension='tsv',
            topic_id=250482
        )
    except IntegrityError:
        db_session.rollback()
        db_session.close()
    	raise ValueError('That file has already been uploaded and cannot be reused. Please change the file contents and try again.')
    tsv_file.seek(0)
    raw_file_content = csv.reader(tsv_file, delimiter='\t', dialect=csv.excel_tab)
    annotations = load_summaries(db_session, raw_file_content, username)
    db_session.close()
    return annotations

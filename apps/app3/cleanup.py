# clean up files in both uploads and results each week.
# This script can be run manually.
# Presently, it is set to run as a cron job each day.
import os
import time

UPLOAD_FOLDER = 'uploads'
RESULT_FOLDER = 'results'
MAX_FILE_AGE = 7 * 24 * 60 * 60  # 7 days in seconds

def cleanup_folder(folder):
    now = time.time()
    for filename in os.listdir(folder):
        filepath = os.path.join(folder, filename)
        if os.path.isfile(filepath):
            file_age = now - os.path.getmtime(filepath)
            if file_age > MAX_FILE_AGE:
                os.remove(filepath)
                print(f"Deleted {filepath}")

if __name__ == '__main__':
    cleanup_folder(UPLOAD_FOLDER)
    cleanup_folder(RESULT_FOLDER)


import sys
from pathlib import Path

if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--self-test':
        from iplan_dicom.smoke import run
        run(Path(sys.argv[2]))
    else:
        from iplan_dicom.gui import main
        main()

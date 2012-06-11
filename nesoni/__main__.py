
import sys

import nesoni
from nesoni import grace, config

try:
    sys.exit(nesoni.main(sys.argv[1:]))

except grace.Help_shown:
    sys.exit(1)

except Exception:
    config.report_exception()
    sys.exit(1) 
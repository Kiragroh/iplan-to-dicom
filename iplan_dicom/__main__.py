"""Command-line entry point; reports omit patient identifiers."""
import argparse
import json
from pathlib import Path
import sys

from .service import inspect_case
from .batch import convert_plan, convert_batch, select_case
from .source import prepare_source
from .xbrain import FormatError


def main(argv=None):
    parser = argparse.ArgumentParser(description='Lokaler iPlan -> DICOM-Konverter')
    commands = parser.add_subparsers(dest='command', required=True)
    inspect = commands.add_parser('inspect', help='Pläne eines Fallordners auflisten')
    inspect.add_argument('case', type=Path)
    inspect.add_argument('--case-index', type=int, help='Fallindex bei ZIP-/Ordnersammlungen (ab 0)')
    convert = commands.add_parser('convert', help='Einen explizit ausgewählten Plan konvertieren')
    convert.add_argument('case', type=Path)
    convert.add_argument('output', type=Path)
    convert.add_argument('--plan-index', type=int, required=True, help='Index aus inspect (ab 0)')
    convert.add_argument('--case-index', type=int, help='Fallindex bei ZIP-/Ordnersammlungen (ab 0)')
    convert.add_argument('--without-mr', action='store_true', help='Nur den CT-basierten RT-Fall ausgeben')
    convert.add_argument('--resume', action='store_true', help='Vorhandene Ausgabe nur nach Quellen- und Hashprüfung übernehmen')
    batch = commands.add_parser('batch', help='Alle aktuellen Pläne eines Falls oder einer Fallsammlung konvertieren')
    batch.add_argument('source', type=Path)
    batch.add_argument('output', type=Path, help='Ausgabewurzel mit Patienten- und Plan-Unterordnern')
    batch.add_argument('--without-mr', action='store_true')
    batch.add_argument('--resume', action='store_true', help='Nur vollständig verifizierte vorhandene Ergebnisse übernehmen')
    commands.add_parser('gui', help='Lokale Oberfläche öffnen')
    args = parser.parse_args(argv)
    try:
        if args.command == 'gui':
            from .gui import main as gui_main
            gui_main()
            return 0
        if args.command == 'inspect':
            with prepare_source(args.case) as source:
                report = inspect_case(select_case(source.root, args.case_index))
        elif args.command == 'batch':
            report = convert_batch(args.source, args.output, include_mr=not args.without_mr, resume=args.resume,
                                   progress=lambda message: print(message, file=sys.stderr, flush=True))
        else:
            report = convert_plan(args.case, args.output, plan_index=args.plan_index,
                                  include_mr=not args.without_mr, resume=args.resume, case_index=args.case_index,
                                  progress=lambda message: print(message, file=sys.stderr, flush=True))
        print(json.dumps(report, ensure_ascii=False, indent=2))
        return 2 if report.get('status') == 'partial_failure' else 0
    except FormatError as exc:
        print(json.dumps({'status': 'stopped', 'reason': str(exc)}, ensure_ascii=False), file=sys.stderr)
        return 2
    except OSError:
        print('Dateizugriff fehlgeschlagen. Quelle und Ausgabeordner prüfen.', file=sys.stderr)
        return 3
    except Exception:
        # Parser/library exception strings may contain source values or paths.
        print(json.dumps({'status': 'stopped', 'reason': 'Unexpected processing error; no complete export confirmed'}),
              file=sys.stderr)
        return 4


if __name__ == '__main__':
    raise SystemExit(main())

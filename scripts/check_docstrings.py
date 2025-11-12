"""
MHKiT-MATLAB Docstring Checker

This script validates MATLAB function docstrings against the MHKiT-MATLAB
docstring template to ensure proper Sphinx documentation generation.
Docstrings are parsed by: https://github.com/sphinx-contrib/matlabdomain
This library works, but the documentation does not provide guidance about
exact docstring formatting requirements. This script checks for common
issues that can cause Sphinx build failures or inconsistent documentation.

This is the expected docstring format for all MHKiT-MATLAB functions:

```
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Brief one-line description of what the function does
%
% Optional extended description paragraph. Can span multiple lines
% to provide additional context, mathematical background, data sources,
% or important behavioral notes about the function.
%
% Parameters
% ------------
% input_1 : type [units]
%   Description of first parameter. Can span multiple lines.
%   Additional details aligned with first line. Use two spaces for indentation
%     input_1.fieldname_1 : type [units]
%       Description of field
%     input_1.fieldname_2 : type [units]
%       Description of field
% input_2 : type [units]
%   Description of second parameter
% input_t : type [units] (optional)
%   Description of second parameter. Mark optional parameters
%   with (optional) after the type.
%
% Returns
% ---------
% output : type
%   Description of return value. For structure outputs,
%   document nested fields with indentation:
%     output.field1 : type [units]
%       Field 1 description
%     output.field2 : type [units]
%       Field description
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
```

This only fails when docstrings are malformed in certain ways that break Sphinx.
Otherwise this shows warnings for incomplete docstrings that should be fixed
but do not break the build.

Exit codes:
    0 - All checks passed (includes warnings)
    1 - Critical errors found (malformed delimiters that break Sphinx builds)

Usage:
    python scripts/check_docstrings.py [--verbose]
"""

import re
import sys
from pathlib import Path
from typing import List, Tuple, Dict
from dataclasses import dataclass
from enum import Enum


class Severity(Enum):
    """Issue severity levels."""

    ERROR = "ERROR"  # Critical - will fail Sphinx build
    WARNING = "WARNING"  # Non-critical - should be fixed but won't break build


@dataclass
class DocstringIssue:
    """Represents a docstring validation issue."""

    file_path: Path
    line_number: int
    severity: Severity
    message: str

    def __str__(self):
        return f"{self.file_path}:{self.line_number}: {self.severity.value}: {self.message}"


class DocstringChecker:
    """Validates MATLAB docstrings against MHKiT standards."""

    # Pattern that caused the Sphinx build failure
    # Only match lines that look like delimiters (4+ consecutive % after initial %)
    MALFORMED_DELIMITER_PATTERNS = [
        re.compile(
            r"^%\s+%{4,}\s*$"
        ),  # % %%%% (space before 4+ percent signs, nothing after)
        re.compile(r"^%%\s+%{4,}\s*$"),  # %% %%%% (space after double percent)
        re.compile(
            r"^%\s*%\s+%{3,}\s*$"
        ),  # % % %%% (spaces between percent signs in delimiter)
    ]

    # Expected patterns
    VALID_DELIMITER = re.compile(r"^%{4,}$")  # 4 or more percent signs, no spaces
    FUNCTION_LINE = re.compile(r"^\s*function\s+")
    DESCRIPTION_LINE = re.compile(
        r"^%\s{1,4}\S+"
    )  # % followed by 1-4 spaces and content
    PARAMETERS_HEADER = re.compile(r"^%\s*Parameters\s*$")
    PARAMETERS_SEP = re.compile(r"^%\s*-{4,}\s*$")
    RETURNS_HEADER = re.compile(r"^%\s*Returns?\s*$")
    RETURNS_SEP = re.compile(r"^%\s*-{4,}\s*$")
    BLANK_COMMENT = re.compile(r"^%\s*$")
    PARAM_LINE = re.compile(r"^%\s{2,}(\w+)\s*:")  # Parameter definition line
    ARGUMENTS_LINE = re.compile(r"^\s*arguments\b")
    END_LINE = re.compile(r"^\s*end\s*$")

    def __init__(self, verbose: bool = False):
        """Initialize the checker.

        Args:
            verbose: If True, print detailed checking information
        """
        self.verbose = verbose
        self.issues: List[DocstringIssue] = []

    def check_file(self, file_path: Path) -> List[DocstringIssue]:
        """Check a single MATLAB file for docstring issues.

        Args:
            file_path: Path to the MATLAB file

        Returns:
            List of issues found in the file
        """
        file_issues = []

        try:
            with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
                lines = f.readlines()
        except Exception as e:
            file_issues.append(
                DocstringIssue(
                    file_path=file_path,
                    line_number=0,
                    severity=Severity.ERROR,
                    message=f"Failed to read file: {e}",
                )
            )
            return file_issues

        # Find function definitions and check their docstrings
        for i, line in enumerate(lines, start=1):
            if self.FUNCTION_LINE.match(line):
                if self.verbose:
                    print(f"  Found function at line {i}")
                issues = self._check_function_docstring(file_path, lines, i - 1)
                file_issues.extend(issues)

        return file_issues

    def _check_function_docstring(
        self, file_path: Path, lines: List[str], func_line_idx: int
    ) -> List[DocstringIssue]:
        """Check the docstring for a function.

        Args:
            file_path: Path to the file being checked
            lines: All lines in the file
            func_line_idx: Index of the function definition line (0-based)

        Returns:
            List of issues found in this function's docstring
        """
        issues = []

        # Skip if this is a nested function (indented)
        if lines[func_line_idx].startswith(" ") or lines[func_line_idx].startswith(
            "\t"
        ):
            if self.verbose:
                print(f"    Skipping nested function at line {func_line_idx + 1}")
            return issues

        # Extract the docstring (from function line to first non-comment line)
        docstring_start = func_line_idx + 1
        docstring_end = docstring_start

        # Skip blank lines after function definition
        while docstring_end < len(lines) and lines[docstring_end].strip() == "":
            docstring_end += 1

        if docstring_end >= len(lines):
            return issues

        # Check if there's a docstring (starts with %)
        if not lines[docstring_end].strip().startswith("%"):
            issues.append(
                DocstringIssue(
                    file_path=file_path,
                    line_number=func_line_idx + 1,
                    severity=Severity.WARNING,
                    message="Function missing docstring",
                )
            )
            return issues

        # Find the end of the docstring
        docstring_start = docstring_end
        while docstring_end < len(lines):
            line = lines[docstring_end].strip()
            if line and not line.startswith("%"):
                break
            docstring_end += 1

        docstring_lines = lines[docstring_start:docstring_end]

        # Remove trailing blank lines from docstring
        while docstring_lines and docstring_lines[-1].strip() == "":
            docstring_lines.pop()

        # Now check the docstring content
        issues.extend(
            self._validate_docstring(file_path, docstring_lines, docstring_start + 1)
        )

        return issues

    def _validate_docstring(
        self, file_path: Path, docstring: List[str], start_line: int
    ) -> List[DocstringIssue]:
        """Validate the content of a docstring.

        Args:
            file_path: Path to the file being checked
            docstring: Lines of the docstring (including % prefix)
            start_line: Line number where docstring starts (1-based)

        Returns:
            List of issues found
        """
        issues = []

        if not docstring:
            return issues

        # Check for malformed delimiters (CRITICAL - causes Sphinx failures)
        for i, line in enumerate(docstring):
            line_num = start_line + i
            for pattern in self.MALFORMED_DELIMITER_PATTERNS:
                if pattern.match(line):
                    issues.append(
                        DocstringIssue(
                            file_path=file_path,
                            line_number=line_num,
                            severity=Severity.ERROR,
                            message=f"Malformed delimiter with space: '{line.rstrip()}' - "
                            f"Must use solid block of % signs (e.g., %%%%%)",
                        )
                    )

        # Check opening delimiter
        first_line = docstring[0].rstrip()
        if not self.VALID_DELIMITER.match(first_line):
            issues.append(
                DocstringIssue(
                    file_path=file_path,
                    line_number=start_line,
                    severity=Severity.WARNING,
                    message=f"Opening delimiter should be 4+ percent signs: '{first_line}'",
                )
            )

        # Check closing delimiter
        last_line = docstring[-1].rstrip()
        if not self.VALID_DELIMITER.match(last_line):
            issues.append(
                DocstringIssue(
                    file_path=file_path,
                    line_number=start_line + len(docstring) - 1,
                    severity=Severity.WARNING,
                    message=f"Closing delimiter should be 4+ percent signs: '{last_line}'",
                )
            )

        # Check for function description (should be after opening delimiter)
        has_description = False
        for i in range(1, min(10, len(docstring))):  # Check first 10 lines
            if self.DESCRIPTION_LINE.match(docstring[i]):
                has_description = True
                break
            if self.PARAMETERS_HEADER.match(docstring[i]):
                break  # Reached parameters without finding description

        if not has_description:
            issues.append(
                DocstringIssue(
                    file_path=file_path,
                    line_number=start_line + 1,
                    severity=Severity.WARNING,
                    message="Missing function description after opening delimiter",
                )
            )

        # Check for Parameters section
        has_parameters = False
        has_parameters_separator = False

        for i, line in enumerate(docstring):
            if self.PARAMETERS_HEADER.match(line):
                has_parameters = True
                # Check if next non-blank line is the separator
                for j in range(i + 1, min(i + 3, len(docstring))):
                    if self.PARAMETERS_SEP.match(docstring[j]):
                        has_parameters_separator = True
                        break
                    if docstring[j].strip() and not self.BLANK_COMMENT.match(
                        docstring[j]
                    ):
                        break
                break

        if not has_parameters:
            issues.append(
                DocstringIssue(
                    file_path=file_path,
                    line_number=start_line,
                    severity=Severity.WARNING,
                    message="Missing 'Parameters' section",
                )
            )
        elif not has_parameters_separator:
            issues.append(
                DocstringIssue(
                    file_path=file_path,
                    line_number=start_line,
                    severity=Severity.WARNING,
                    message="Missing separator (----) under 'Parameters' header",
                )
            )

        # Check for Returns section
        has_returns = False
        has_returns_separator = False

        for i, line in enumerate(docstring):
            if self.RETURNS_HEADER.match(line):
                has_returns = True
                # Check if next non-blank line is the separator
                for j in range(i + 1, min(i + 3, len(docstring))):
                    if self.RETURNS_SEP.match(docstring[j]):
                        has_returns_separator = True
                        break
                    if docstring[j].strip() and not self.BLANK_COMMENT.match(
                        docstring[j]
                    ):
                        break
                break

        if not has_returns:
            issues.append(
                DocstringIssue(
                    file_path=file_path,
                    line_number=start_line,
                    severity=Severity.WARNING,
                    message="Missing 'Returns' section",
                )
            )
        elif not has_returns_separator:
            issues.append(
                DocstringIssue(
                    file_path=file_path,
                    line_number=start_line,
                    severity=Severity.WARNING,
                    message="Missing separator (----) under 'Returns' header",
                )
            )

        return issues

    def check_directory(self, directory: Path) -> None:
        """Check all MATLAB files in a directory recursively.

        Args:
            directory: Directory to check
        """
        matlab_files = sorted(directory.rglob("*.m"))

        if self.verbose:
            print(f"Checking {len(matlab_files)} MATLAB files in {directory}...")

        for file_path in matlab_files:
            if self.verbose:
                print(f"\nChecking {file_path}...")
            issues = self.check_file(file_path)
            self.issues.extend(issues)

    def print_summary(self) -> Tuple[int, int]:
        """Print a summary of all issues found.

        Returns:
            Tuple of (error_count, warning_count)
        """
        errors = [issue for issue in self.issues if issue.severity == Severity.ERROR]
        warnings = [
            issue for issue in self.issues if issue.severity == Severity.WARNING
        ]

        # Group by file for better readability
        issues_by_file: Dict[Path, List[DocstringIssue]] = {}
        for issue in self.issues:
            if issue.file_path not in issues_by_file:
                issues_by_file[issue.file_path] = []
            issues_by_file[issue.file_path].append(issue)

        if self.issues:
            print("\n" + "=" * 80)
            print("DOCSTRING CHECK RESULTS")
            print("=" * 80 + "\n")

            for file_path in sorted(issues_by_file.keys()):
                print(f"\n{file_path}:")
                for issue in issues_by_file[file_path]:
                    print(
                        f"  Line {issue.line_number}: [{issue.severity.value}] {issue.message}"
                    )

        print("\n" + "=" * 80)
        print(f"Total: {len(errors)} errors, {len(warnings)} warnings")
        print("=" * 80)

        if errors:
            print("\nFAILED: Critical errors found that will break Sphinx build!")
            print("   Fix all ERROR issues before committing.")
        elif warnings:
            print("\nWARNINGS: Some docstrings are incomplete or inconsistent.")
            print("   Please review and fix WARNING issues.")
        else:
            print("\nPASSED: All docstrings are properly formatted!")

        return len(errors), len(warnings)


def main():
    """Main entry point for the script."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Check MHKiT-MATLAB docstrings for format compliance"
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Print verbose checking information",
    )
    parser.add_argument(
        "--path",
        type=Path,
        default=Path("mhkit"),
        help="Path to check (default: mhkit)",
    )

    args = parser.parse_args()

    # Resolve path relative to script location
    script_dir = Path(__file__).parent.parent
    check_path = script_dir / args.path

    if not check_path.exists():
        print(f"Error: Path does not exist: {check_path}", file=sys.stderr)
        return 1

    print(f"Checking docstrings in: {check_path}")

    checker = DocstringChecker(verbose=args.verbose)
    checker.check_directory(check_path)
    error_count, warning_count = checker.print_summary()

    # Determine exit code
    # Exit 1 only for critical errors that break Sphinx builds
    # Warnings don't fail the build
    if error_count > 0:
        return 1
    else:
        return 0


if __name__ == "__main__":
    sys.exit(main())

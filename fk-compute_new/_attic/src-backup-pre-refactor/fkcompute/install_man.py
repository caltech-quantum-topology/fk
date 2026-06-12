#!/usr/bin/env python3
"""
Install the fk man page to the system.

This script attempts to install the fk.1 man page to the appropriate
system location and update the man database.
"""

import os
import sys
import shutil
import subprocess
import platform
from pathlib import Path
from importlib import resources


def find_man_page() -> Path:
    """Find the man page file in the package or source directory."""
    try:
        # Try to find it in the installed package resources
        with resources.path('fkcompute', '__init__.py') as pkg_path:
            man_file = pkg_path.parent.parent.parent / "fk.1"
            if man_file.exists():
                return man_file
    except Exception:
        pass

    # Try to find it relative to this script (development)
    script_dir = Path(__file__).parent
    for path in [
        script_dir / "../../fk.1",          # From src/fkcompute/
        script_dir / "../../../fk.1",       # Alternative layout
        Path.cwd() / "fk.1",               # Current directory
    ]:
        if path.exists():
            return path.resolve()

    raise FileNotFoundError("Could not find fk.1 man page file")


def get_man_directories():
    """Get list of potential man page directories in order of preference."""
    directories = []

    if platform.system() in ('Linux', 'Darwin'):  # Unix-like systems
        # User-local first (no sudo required)
        home_man = Path.home() / ".local/share/man/man1"
        directories.append(home_man)

        # System-wide directories
        directories.extend([
            Path("/usr/local/share/man/man1"),
            Path("/usr/share/man/man1"),
        ])

        # Check for Homebrew on macOS
        if platform.system() == 'Darwin':
            brew_prefix = os.environ.get('HOMEBREW_PREFIX')
            if not brew_prefix:
                try:
                    result = subprocess.run(['brew', '--prefix'],
                                          capture_output=True, text=True, timeout=5)
                    if result.returncode == 0:
                        brew_prefix = result.stdout.strip()
                except (subprocess.SubprocessError, FileNotFoundError):
                    pass

            if brew_prefix:
                directories.insert(-1, Path(brew_prefix) / "share/man/man1")

    return directories


def install_man_page(source: Path, destination: Path) -> bool:
    """Install man page to destination directory."""
    try:
        # Create directory if it doesn't exist
        destination.mkdir(parents=True, exist_ok=True)

        # Test write permission
        test_file = destination / ".test_write"
        try:
            test_file.touch()
            test_file.unlink()
        except (PermissionError, OSError):
            return False

        # Copy the man page
        man_dest = destination / "fk.1"
        shutil.copy2(source, man_dest)

        print(f"✓ Man page installed to: {man_dest}")
        return True

    except (PermissionError, OSError, Exception) as e:
        print(f"✗ Failed to install to {destination}: {e}")
        return False


def update_man_database(man_dir: Path) -> None:
    """Update the man database."""
    try:
        # Try to update man database
        commands = ['mandb', 'makewhatis']  # Different systems use different commands

        for cmd in commands:
            try:
                # Try user-local update first
                if str(man_dir).startswith(str(Path.home())):
                    subprocess.run([cmd, str(man_dir.parent)],
                                 check=False, stdout=subprocess.DEVNULL,
                                 stderr=subprocess.DEVNULL, timeout=30)
                else:
                    # System-wide update
                    subprocess.run([cmd], check=False,
                                 stdout=subprocess.DEVNULL,
                                 stderr=subprocess.DEVNULL, timeout=30)
                print(f"✓ Man database updated with {cmd}")
                return
            except (subprocess.SubprocessError, FileNotFoundError):
                continue

        print("! Note: Could not update man database automatically")
        print("  Run 'mandb' or 'makewhatis' to update manually")

    except Exception as e:
        print(f"! Warning: Failed to update man database: {e}")


def main():
    """Main installation function."""
    print("Installing fk man page...")

    # Check if we're on a supported system
    if platform.system() not in ('Linux', 'Darwin'):
        print(f"✗ Man page installation not supported on {platform.system()}")
        print("  Man pages are primarily for Unix-like systems")
        return 1

    try:
        # Find the man page file
        man_source = find_man_page()
        print(f"Found man page: {man_source}")

        # Try installation directories in order
        directories = get_man_directories()
        installed = False

        for man_dir in directories:
            print(f"Trying: {man_dir}")
            if install_man_page(man_source, man_dir):
                update_man_database(man_dir)
                installed = True
                break

        if not installed:
            print("✗ Failed to install man page to any directory")
            print("\nManual installation:")
            print(f"  sudo cp {man_source} /usr/local/share/man/man1/")
            print("  sudo mandb")
            return 1

        print("\n✓ Man page installation complete!")
        print("  You can now use 'man fk' to view the manual")
        return 0

    except FileNotFoundError as e:
        print(f"✗ Error: {e}")
        return 1
    except Exception as e:
        print(f"✗ Unexpected error: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
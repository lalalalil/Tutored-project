#!/usr/bin/env python3
"""
Cleanup Script: Keep ONLY SBML models + JSON metadata
======================================================

Removes all files except .sbml and .json from the downloaded articles.
Keeps folder structure intact.

Usage:
    python cleanup_files.py               # Dry run (shows what will be deleted)
    python cleanup_files.py --delete      # Actually delete files
"""

import logging
import shutil
from pathlib import Path
from collections import defaultdict

logging.basicConfig(level=logging.INFO, format="%(levelname)s | %(message)s")
logger = logging.getLogger("cleanup")

# KEEP ONLY THESE EXTENSIONS
KEEP_EXTENSIONS = {".sbml", ".json"}

def get_output_root() -> Path:
    """Get downloaded_ncbi_models folder."""
    home = Path.home()
    return home / "Desktop" / "downloaded_ncbi_models"


def cleanup_directory(root_dir: Path, dry_run: bool = False) -> dict:
    """
    Delete all files except .sbml and .json.

    Returns stats dict with file counts.
    """
    stats = {
        "files_deleted": 0,
        "files_kept": 0,
        "total_size_deleted": 0,
        "by_extension": defaultdict(int),
        "errors": 0,
    }

    logger.info("Scanning: %s", root_dir)

    # Find all files
    for path in sorted(root_dir.rglob("*")):
        if not path.is_file():
            continue

        suffix = path.suffix.lower()

        if suffix in KEEP_EXTENSIONS:
            # KEEP THIS FILE
            stats["files_kept"] += 1
            logger.info("  ✓ KEEP: %s", path.name)

        else:
            # DELETE THIS FILE
            size = path.stat().st_size
            stats["files_deleted"] += 1
            stats["total_size_deleted"] += size
            stats["by_extension"][suffix] += 1

            if dry_run:
                logger.info("  🗑️  DELETE: %s (%d bytes)", path.name, size)
            else:
                try:
                    path.unlink()
                    logger.info("  🗑️  DELETED: %s (%d bytes)", path.name, size)
                except Exception as e:
                    stats["errors"] += 1
                    logger.warning("  ❌ ERROR: %s", e)

    # Delete empty directories
    for path in sorted(root_dir.rglob("*"), reverse=True):
        if path.is_dir() and not any(path.iterdir()):
            if dry_run:
                logger.info("  📁 DELETE (empty): %s", path.name)
            else:
                try:
                    path.rmdir()
                    logger.info("  📁 DELETED (empty): %s", path.name)
                except Exception as e:
                    logger.warning("  ❌ ERROR: %s", e)

    return stats


def cleanup_all_diseases(dry_run: bool = False):
    """Clean up all disease directories."""
    output_root = get_output_root()

    if not output_root.exists():
        logger.error("Not found: %s", output_root)
        return

    logger.info("\n" + "="*70)
    logger.info("CLEANUP: %s", output_root)
    logger.info("="*70)

    total_stats = {
        "files_deleted": 0,
        "files_kept": 0,
        "total_size_deleted": 0,
        "by_extension": defaultdict(int),
        "errors": 0,
    }

    # Clean each disease folder
    for disease_dir in sorted(output_root.iterdir()):
        if not disease_dir.is_dir() or disease_dir.name.startswith("."):
            continue

        logger.info("\n📁 %s", disease_dir.name)

        # Clean each PMC article
        for pmc_dir in disease_dir.iterdir():
            if not pmc_dir.is_dir() or pmc_dir.name.startswith("_"):
                continue

            # Clean data folder
            data_dir = pmc_dir / "data"
            if data_dir.exists():
                stats = cleanup_directory(data_dir, dry_run=dry_run)
                total_stats["files_deleted"] += stats["files_deleted"]
                total_stats["files_kept"] += stats["files_kept"]
                total_stats["total_size_deleted"] += stats["total_size_deleted"]
                total_stats["errors"] += stats["errors"]
                for ext, count in stats["by_extension"].items():
                    total_stats["by_extension"][ext] += count

            # Clean metadata folder
            metadata_dir = pmc_dir / "metadata"
            if metadata_dir.exists():
                stats = cleanup_directory(metadata_dir, dry_run=dry_run)
                total_stats["files_deleted"] += stats["files_deleted"]
                total_stats["files_kept"] += stats["files_kept"]
                total_stats["total_size_deleted"] += stats["total_size_deleted"]
                total_stats["errors"] += stats["errors"]
                for ext, count in stats["by_extension"].items():
                    total_stats["by_extension"][ext] += count

    # Print summary
    logger.info("\n" + "="*70)
    logger.info("SUMMARY")
    logger.info("="*70)

    logger.info("\n✓ Files KEPT: %d", total_stats['files_kept'])
    logger.info("  Files DELETED: %d", total_stats['files_deleted'])
    logger.info(" Space freed: %.2f MB", total_stats['total_size_deleted'] / (1024**2))
    logger.info(" Errors: %d", total_stats['errors'])

    if total_stats["by_extension"]:
        logger.info("\nDeleted by type (top 10):")
        for ext, count in sorted(total_stats["by_extension"].items(), key=lambda x: -x[1])[:10]:
            logger.info("  %s: %d files", ext, count)

    if dry_run:
        logger.info("\n  DRY RUN: Nothing was deleted!")
        logger.info("   Run with --delete to actually delete:\n")
        logger.info("   python cleanup_files.py --delete\n")
    else:
        logger.info("\n Cleanup complete!")


if __name__ == "__main__":
    import sys

    dry_run = "--delete" not in sys.argv

    if dry_run:
        logger.info("\n" + "="*70)
        logger.info("  DRY RUN MODE (no files will be deleted)")
        logger.info("="*70 + "\n")

    cleanup_all_diseases(dry_run=dry_run)

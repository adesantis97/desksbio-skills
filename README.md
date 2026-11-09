# desks.bio skills

Free, open-source agent skills for biologists, built at the bench and used on real experiments.
Each folder is one skill in the open `SKILL.md` format: it works in Claude Code, Cowork and any agent that reads skills.

Every skill keeps a human in the loop: the agent does the tedious part, you check it on a page and correct with a click.

## Skills

- **colony-count**: colony formation assay from photos of the dishes. Colonies per dish, mean colony size, fold change vs CTRL.
- **invasion-count**: Transwell invasion/migration from microscope fields. Invaded cells per field (one dot per cell, click to fix), invaded area beside, fold change across biological replicates, Excel table.

## Install

Copy one skill into your skills folder:

```bash
git clone https://github.com/adesantis97/desksbio-skills.git
cp -r desksbio-skills/colony-count ~/.claude/skills/
```

Python 3.10 to 3.13, once:

```bash
pip install -r desksbio-skills/requirements.txt
```

R + tidyverse are needed only if you ask for figures.

Then drop your photos in a folder and say *"count the colonies"* or *"measure this invasion assay"*.

Prefer a Claude Code plugin? The same skills ship as `bench-assays` in [adesantis97/desks](https://github.com/adesantis97/desks).

## More

New skills come one at a time, each with a short video on [@desks.bio](https://instagram.com/desks.bio).
Built by the team behind LabsDesk and OmicsDesk.

## License

MIT. See [LICENSE](LICENSE).

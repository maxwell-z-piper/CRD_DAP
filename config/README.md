# Target configuration files

Copy `target_config_template.py` to one file per target and edit that copy.

Example:

```bash
cp config/target_config_template.py config/9028-1901.py
```

The target file contains both required input paths and configurable numerical-analysis settings. Development defaults are labeled clearly and should not be treated as universal scientific values.

Every pipeline run saves a verbatim `config_snapshot.py` and JSON configuration manifest inside the run metadata directory.

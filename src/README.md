# src

Python implementations for the WPMW project.

## Output path convention

All scripts in this directory must write files through `output_path()` from
`wpmw_utils`, never via hardcoded absolute paths:

```python
from wpmw_utils import output_path

fig.savefig(output_path("my_figure.png"), dpi=150, bbox_inches="tight")
```

The output directory is controlled by the `WPMW_OUTPUT` environment variable.
If unset, files go to `./output`.

See the top-level `README.md` for the full convention.

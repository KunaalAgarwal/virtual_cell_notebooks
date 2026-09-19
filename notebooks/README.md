# notebooks/

Exploration and visualization only.

Model logic never lives here, even temporarily — it goes in `models/` as a `.py` file and gets imported:

```python
from models.baseline_control_mean import predict
```

Existing research notebooks (published-approach reimplementations, data pulls) remain in `src/`.

# ML Design Pattern

- **Date:** 30-Mar-2021
- **Presenter:** Valliappa Lakshmanan
- **Event:** Meetup - Big Data Madison

Github Page for Book: **Machine Learning Design Patterns** - [github](https://github.com/GoogleCloudPlatform/ml-design-patterns "Github @GoogleCloudPlatform/ml-design-patterns")

- Book: Machine Learning Design Patterns
- Authors:
  - Valliappa Lakshmanan
  - Sara Robinson
  - Michael Munn

## 02 Add Neutral Class

Employ neutral class to:

    - improve embedding with confident classifications

## 03 Checkpoints

    1. Resilience during long training times
    2. Generalization (early stopping)
        - Better approach is to add regularization
    3. Fine-tuning
        - model trained -> production -> new data -> start training from earlier checkpoint
            e.g.

```python
# Using Keras
ckpt_callback = tf.keras.callbacks.ModelCheckpoint(...)
history = model.fit(..., callback=[ckpt_callback])
```

Redefine what an epoch is (i.e. _NOT_ an integer epoch)

```python
trainds = trainds.repeat()
# Set num of epochs to be the num of checkpoints
NUM_STEPS = 143000
BACTCH_SIZE = 50
NUM_CHECKPOINTS = 15
cp_callback = ...
history = model.fit(trainds,
validation_data=evalds,
epocs=NUM_CHECKPOINTS,
steps_per_epoch=steps_per_epoch,
batch_size=BATCH_SIZE
)
```

# 04 Pass-through keys in Keras

- Keyed predictions --> batch predictions & async predictions
- Sliced evaluation (How are we slicing our training data?)

# 05 Transform

Transform ensures transformations are automatically applied during `ML.PREDICT`

In TF/Keras, do transformations in Lambda Layers so that they are part of the model graph.

Transform pattern: The model graph should include the transformations

Do "Feature Store" if pre-processing is becoming majorly complex

# 06 Summary

1. Feature Cross
2. Neutral Class
3.

Machine Learning Design Patterns

https://bit.ly/ml-design-patterns

## 07 Extra

Can you treat model as network of models? Yes, but how do you track dependencies and artifacts?

Cascade -- Model Pipeline

## Model and Data Versioning best practices

Config file to configure each container

[Kubeflow Pipelines](kubeflow.org) for ML model data versioning

> ### Note:
>
> Troy:
>
> feature cross really only makes sense to apply to the independent / predictor variables, or perhaps among dependent / response variables if you have multiple of those. Any other crosses would mix inputs and outputs.

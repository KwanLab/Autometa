from django.db import models
from django.utils import timezone
from django.contrib.auth.models import User
# Object relational mapper (ORM) allows access to database in an OOP way. You
# can change databases (like SQLite vs Postgres) by setting it up in the
# settings. The queries will remain the same. We can represent our database
# structure using classes but here in django, they are called models.


class Project(models.Model):
    title = models.CharField(max_length=127)
    date_created = models.DateTimeField(default=timezone.now)
    description = models.TextField()
    # link using many-to-one relationship to a user (using the built-in User class)
    user = models.ForeignKey(User, on_delete=models.CASCADE)

    # adding the Meta changes how the tables are created as explained here:
    # https://stackoverflow.com/a/16655437 Django will create tables from the
    # classes that inherit from Project regardless of the fields defined
    # (attributes).
    # class Meta:
    #     abstract = True

    def __str__(self):
        return self.title


class Job(models.Model):
    title = models.CharField(max_length=127)
    date_run = models.DateTimeField(default=timezone.now)
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    # link a Job to a Project using a many-to-one relationship
    project = models.ForeignKey(
        Project, on_delete=models.CASCADE, null=True, blank=False)

    def __str__(self):
        return self.title

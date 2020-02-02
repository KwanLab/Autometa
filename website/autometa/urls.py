from django.urls import path
from . import views

urlpatterns = [
    path('', views.startpage, name='website-startpage'),
]

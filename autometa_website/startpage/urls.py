from django.urls import path
from . import views

urlpatterns = [
    path('', views.startpage, name='startpage'),
    path('about/', views.about, name='about'),
]

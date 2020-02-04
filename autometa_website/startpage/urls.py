from django.urls import path
from . import views

urlpatterns = [
    path('', views.startpage, name='startpage-home'),
    path('about/', views.about, name='startpage-about'),
]

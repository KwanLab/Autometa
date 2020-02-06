from django.shortcuts import render
from .models import Project


def startpage(request):
    context = {'project': Project.objects.all()}
    return render(request, 'startpage/home.html', context)


def about(request):
    return render(request, 'startpage/about.html', context={'title': 'about autometa'})

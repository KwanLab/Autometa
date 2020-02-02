from django.shortcuts import render
from django.http import HttpResponse

# Create your views here.


def startpage(request):
    return HttpResponse('<h1>Welcome to Autometa</h1>')

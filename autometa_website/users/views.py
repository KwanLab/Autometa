from django.shortcuts import render, redirect
from django.contrib import messages
from django.contrib.auth.decorators import login_required
from .forms import UserRegisterForm

# Create your views here.


def register(request):
    # If [if] a form is submitted by a user, i.e. a POST request is received, we take
    # in the request, validate it and save it. Otherwise [else], if the form is
    # requested, then we send back a blank form
    if request.method == 'POST':
        form = UserRegisterForm(request.POST)
        if form.is_valid():
            # saves the user
            form.save()
            # valid form data will be in the dictionary called 'form.cleaned_data'
            username = form.cleaned_data.get('username')
            # show flash message to show that we've received the data
            messages.success(
                request, f'Your account has been created! You are now able to log in.')
            # then redirect the user to the home page
            return redirect('login')
    else:
        form = UserRegisterForm()
    return render(request, 'users/register.html', {'form': form})


def profile(request):
    return render(request, 'users/profile.html')

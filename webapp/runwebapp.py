from app import app

if __name__ == '__main__':
    # Turn off debug for production environment
    app.run(debug=True,passthrough_errors=False)
    #app.run(debug=False)

import random
import string

def generate_password(length: int) -> str:
    chars = string.ascii_letters + string.digits + string.punctuation
    return ''.join(random.choice(chars) for i in range(length))


password = generate_password(12)
print(password)






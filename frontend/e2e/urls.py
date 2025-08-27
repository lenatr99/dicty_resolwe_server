base_urls = {
    "dev": "http://localhost:3000/",
    "qa": "https://qa.dictyexpress.org/",
    "app": "https://app.dictyexpress.org/",
}


def get_urls(env):
    return {
        "base": base_urls[env],
        "root": f"{base_urls[env]}bcm/",
        "links": {
            "genialis": "https://www.genialis.com/",
            "terms_of_service": "https://www.genialis.com/terms-of-service/",
            "privacy_policy": "https://www.genialis.com/privacy-policy/",
            "help": "https://genialis.zendesk.com/hc/en-us",
        },
    }
